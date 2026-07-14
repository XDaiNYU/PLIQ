# PLIQ ver5 — 各 term 详细计算过程

本文档记录精简后 PLIQ ver5 实际产出的全部分数：**eDockQ1** 与 **eDockQ6_1 / 6_2 / 6_3 / 6_4**，以及每个 term 中 OTMol、原子、残基的映射方式。

---

## 0. 唯一的配体映射：OTMol（只算一次，全程复用）

每个 (case, model) 只调用一次 `align_ligands_otmol_bci_constrained(mol_ref, mol_dock)`，得到一张映射表：

```
dock_to_ref :  docked 配体 RDKit 原子 index  →  reference 配体 RDKit 原子 index
```

这张表在 `runner.py` 计算一次后，被下列所有用到配体原子的环节复用：
- 配体 RMSD（k_l）
- 界面 RMSD（k_i）的配体部分
- FNAT（fnat_BBSC）的配体部分
- BINANA 四个 level 的配体签名部分

**没有 MCS**。蛋白侧不经过 OTMol：ref 与 dock 蛋白共享同一套残基编号/残基名/原子名，直接按氨基酸身份字符串匹配。

辅助映射 `binana_pdbindex_to_rdkit`（BINANA 的 PDB 序号 → dock RDKit index，按 3D 坐标 + 元素就近一一配对）只是把 BINANA 输出的原子接回到 OTMol 的 index 空间，不是第二次比对。

---

## 1. eDockQ1 = (fnat_BBSC + k_i + k_l) / 3

三个子项，三种映射。

### 1a. fnat_BBSC（天然接触保留率，BBSC F1）
函数 `calculate_fnat_backbone_sidechain_otmol`。

- **配体原子映射**：遍历 `dock_to_ref` 每一对 `(dock_idx, ref_idx)`（OTMol）。
- **蛋白/残基映射**：对每个配体原子，在各自结构内找 5 Å 内的蛋白原子，生成接触标签
  `mol_ref{ref_idx}_res{resnum}_{resname}[_backbone/_sidechain]`。
  - 配体部分永远用 `ref_idx`（OTMol 后的天然原子号）→ ref/dock 两集合可直接取交集。
  - 蛋白部分按残基号 + 残基名标识（无需 index 映射）；按原子名分 backbone(N/CA/C/O) 与 sidechain。
- **BBSC**：backbone 接触集 ∪ sidechain 接触集，`F1 = 2·TP / (2·TP + FP + FN)`。

### 1b. k_i = 1 / (1 + (iRMSD / 1.5)²)（界面 RMSD 核）
函数 `calculate_interface_rmsd_otmol`，**不做叠合**。

- **蛋白侧映射**：两侧各取"任一原子距配体 10 Å 内"的界面残基，取交集；按 `(resnum, resname, 原子名)` 配对（默认 CA）。
- **配体侧映射**：界面配体原子在天然侧定义（距 ref 蛋白 10 Å 内），用 `dock_to_ref` 反查 dock 对应原子。
- 合并蛋白 + 配体配对坐标算 RMSD。

### 1c. k_l = 1 / (1 + (ligRMSD / 8.5)²)（配体 RMSD 核）
直接对 `dock_to_ref` 全部配体原子对求 RMSD（无叠合，纯 OTMol）。

> 示例 1B9J model 0：fnat=0.8791，iRMSD=0.5392→k_i=0.8856，ligRMSD=0.9306→k_l=0.9882，**eDockQ1=0.9176**。

---

## 2. eDockQ6_k = term_4 × (B6·H + k_i + k_l)/3 × TM   (k = 1..4)

`k_i`、`k_l` 与 eDockQ1 相同。新增 B6、H、term_4、TM。

### 2a. B6 = 六类 BINANA 合并成一项的 micro-F1（**不是连乘**）
六类非疏水相互作用：`hbond, salt, pipi, tstack, cationpi, halogen` 被**当作一项**处理：先把六类的 TP、ref_n、dock_n **分别求和**，再用合并后的计数算唯一的 precision / recall / F1。
- precision = Σtp / Σdock_n，recall = Σtp / Σref_n，F1 = 2PR/(P+R)。
- 六类在 ref 中全都不存在（Σref_n=0）→ P=R=F1=1.0（无可复现，记为中性 1）。
- Σref_n>0 但 Σdock_n=0，或 Σtp=0 → F1=0。
- 注意：这是 micro 池化，**不会**因为某一类 F1=0 就把整项清零；只有当六类合并后一个 TP 都没有时才为 0。

### 2b. H = hydrophobic F1（单独一项，乘上去）
`B6·H` = 六类合并 micro-F1 × 疏水 F1。

### 2c. BINANA 签名与匹配
每条相互作用编码为 `L[配体部分]R[受体部分]`，ref multiset 与 dock multiset 取逐项 min 交集得 TP。

**配体部分（四个 level 完全相同，均经 OTMol）**：
- ref 侧：BINANA 原子 →（坐标 map `binana_pdbindex_to_rdkit_ref`）→ `mol_ref` 的天然 RDKit index。
- dock 侧：BINANA 原子 → dock RDKit index →（OTMol `dock_to_ref`）→ 天然 RDKit index。
- 即：**dock 原子只要被 OTMol 对到 ref 的对应配体原子即视为同一原子**。签名里写成天然 RDKit 整数 index（如 `L[7|8]`）。
- **氢原子忽略**：RDKit 配体分子按 `removeHs=True` 载入，氢没有 RDKit/OTMol 对应物，因此签名只保留重原子（H-键/盐桥只按重原子供体/受体配对）。
- **ref 侧坐标 map（关键修复）**：BINANA 会把 HETATM 配体改标成 chain `X`/`UNK`，而 RDKit `mol_ref` 是 chain `B`，纯靠 (chain,resID,resName,atomName) 元组匹配会失败 → 旧版 ref 签名静默回退到原始 PDB 序号，永远配不上 dock 的 RDKit 整数签名，导致小分子配体 BINANA 召回几乎全为 0（肽类也被低估）。现在 ref 与 dock 一样用坐标 map 投到 RDKit index，两侧落在同一索引空间。

**受体部分（四个 level 唯一的区别，不经 OTMol）**：

| level | 列 | 受体标签 | 粒度 |
|---|---|---|---|
| eDockQ6_1 | default | chain : resID : resName : atomName | 原子级 |
| eDockQ6_2 | rxatm | chain : resID : resName : atomName | 原子级（当前与 6_1 等价）|
| eDockQ6_3 | rxidx | chain : resID : resName(统一大写) : atomName | 原子级，大小写不敏感 |
| eDockQ6_4 | rxres | chain : resID : resName（丢 atomName）| 残基级，最宽松 |

**level 差异示例（1B9J model 45，hydrophobic）**：同一条配体原子 15 ↔ GLU32 接触，ref 碰 CB/CG、dock 碰 CD：
- 原子级 6_1/6_2/6_3：`L[15]R[A:32:GLU:CB]` ≠ `L[15]R[A:32:GLU:CD]` → 不匹配。
- 残基级 6_4：`L[15]R[A:32:GLU]` 两侧相同 → 匹配。
结果该 level F1 从 0.324 升到 0.378。

### 2d. term_4 = (PoseBusters 8 项检查通过比例)²
8 项：bond_lengths、bond_angles、aromatic_ring_flatness、double_bond_flatness、internal_steric_clash、energy_ratio≤100、minimum_distance_to_protein>0.75、volume_overlap_with_protein。

### 2e. TM = 蛋白 TMscore

> 示例 1B9J model 0（term_4=1.0，TM=1.0 为占位值）：六类合并计数 Σtp=6、Σref=6、Σdock=7 → B6 micro-F1=12/13=0.9231（hbond F1=0.909、salt F1=1.0）；H(hydrophobic)=0.4865 → B6·H=0.4491。四个 level 同为 `1.0 × (0.4491 + 0.8856 + 0.9882)/3 × 1.0 = 0.7743`。
>
> 示例 T1124 model 0（小分子配体，term_4/TM 为占位值）：ligand RMSD 0.45 但 interface RMSD 2.28，配体相对蛋白被整体平移 → 多数原生相互作用未复现。hbond 1/1/7/1（recall=1 但 6 个 FP），hydrophobic F1=0.1587，B6 micro-F1=0.1818 → B6·H=0.0289，eDockQ6 ≈ `1.0 × (0.0289 + 0.302 + 0.9972)/3 = 0.4427`。

---

## 3. 实现要点（ver5 精简后）

- 配体映射 OTMol **只算一次**（`runner.py` 每个 (case, model)），下游全部复用同一张 `dock_to_ref`。
- **无 MCS**：代码中已无 `FindMCS / rdFMCS / GetSubstructMatch`，仅保留"no MCS"说明性注释。
- **不静默**：依赖缺失/计算失败的 case/model 不会被伪造成功，而是记录到 `e-DockQ_results_failures.csv` / `posebusters_results_failures.csv`（含失败原因）。
- 打分函数 `score_edockq_stepwise` 只产出 eDockQ1 与 eDockQ6_1..4，已删除 eDockQ4/5/7–15 等未使用分支；B6 采用六类 micro 池化（`add_binana_six_pool_columns`），不再连乘。

## 4. 复现脚本

- `scripts/run_one_case_terms.py` — 单 case 打印所有中间项与最终分数。
- `scripts/show_level_signatures.py [CASE MODEL]` — 打印四个 level 的 BINANA 签名与汇总表，直观显示 level 差异。
