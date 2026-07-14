# PLIQ ver5 — Detailed per-term calculation

This documents every score the trimmed PLIQ ver5 produces: **eDockQ1** and **eDockQ6_1 / 6_2 / 6_3 / 6_4**, plus how OTMol, atoms, and residues are mapped in each term.

---

## 0. The single ligand mapping: OTMol (computed once, reused everywhere)

For each (case, model), `align_ligands_otmol_bci_constrained(mol_ref, mol_dock)` is called **once**, producing one table:

```
dock_to_ref :  docked ligand RDKit atom index  →  reference ligand RDKit atom index
```

Computed once in `runner.py`, this table is reused by everything that touches ligand atoms:
- ligand RMSD (k_l)
- the ligand part of interface RMSD (k_i)
- the ligand part of FNAT (fnat_BBSC)
- the ligand part of every BINANA signature (all four levels)

**No MCS** anywhere. The protein side does not use OTMol: reference and docked proteins share the same residue numbering/name/atom name, so they are matched by amino-acid identity strings.

The auxiliary map `binana_pdbindex_to_rdkit` (BINANA PDB serial → dock RDKit index, nearest-by-coordinate + element, one-to-one) only re-attaches BINANA's atoms into the OTMol index space; it is not a second alignment.

---

## 1. eDockQ1 = (fnat_BBSC + k_i + k_l) / 3

### 1a. fnat_BBSC (fraction of native contacts, BBSC F1)
`calculate_fnat_backbone_sidechain_otmol`.
- **Ligand mapping**: iterate each `(dock_idx, ref_idx)` in `dock_to_ref` (OTMol).
- **Residue/protein mapping**: for each ligand atom, find protein atoms within 5 Å in each structure; build a contact tag `mol_ref{ref_idx}_res{resnum}_{resname}[_backbone/_sidechain]`. The ligand part always uses `ref_idx` so ref/dock sets intersect directly. Protein side keyed by residue number + name; split backbone (N/CA/C/O) vs sidechain by atom name.
- **BBSC**: backbone-set ∪ sidechain-set, `F1 = 2·TP/(2·TP+FP+FN)`.

### 1b. k_i = 1/(1+(iRMSD/1.5)²) (interface RMSD kernel)
`calculate_interface_rmsd_otmol`, **no superposition**.
- **Protein side**: interface residues (any atom within 10 Å of ligand) on both sides, intersected; paired by `(resnum, resname, atomName)` (CA by default).
- **Ligand side**: interface ligand atoms defined on the reference side, mapped to dock via `dock_to_ref`.

### 1c. k_l = 1/(1+(ligRMSD/8.5)²) (ligand RMSD kernel)
RMSD over all `dock_to_ref` atom pairs (no superposition; OTMol only).

> Example 1B9J model 0: fnat=0.8791, iRMSD=0.5392→k_i=0.8856, ligRMSD=0.9306→k_l=0.9882, **eDockQ1=0.9176**.

---

## 2. eDockQ6_k = term_4 × (B6·H + k_i + k_l)/3 × TM   (k = 1..4)

`k_i`, `k_l` identical to eDockQ1.

### 2a. B6 = ONE micro-pooled F1 over the six BINANA types (**not a product**)
The six non-hydrophobic types `hbond, salt, pipi, tstack, cationpi, halogen` are treated as **a single term**: sum their TP, ref_n and dock_n separately, then compute one precision / recall / F1 from the pooled counts.
- precision = Σtp / Σdock_n, recall = Σtp / Σref_n, F1 = 2PR/(P+R).
- None of the six present in ref (Σref_n=0) → P=R=F1=1.0 (nothing to reproduce; neutral 1).
- Σref_n>0 with Σdock_n=0, or Σtp=0 → F1=0.
- Because this is micro pooling, a single type with F1=0 does **not** zero the whole term; B6=0 only when there is not a single TP across all six pooled.

### 2b. H = hydrophobic F1 (separate factor)
`B6·H` = six-type pooled micro-F1 × hydrophobic F1.

### 2c. BINANA signatures and matching
Each interaction → `L[ligand]R[receptor]`; TP = per-element min intersection of the ref and dock multisets.

**Ligand part (identical across all four levels, all via OTMol)**:
- ref: BINANA atom → (coordinate map `binana_pdbindex_to_rdkit_ref`) → reference RDKit index on `mol_ref`.
- dock: BINANA atom → dock RDKit index → (OTMol `dock_to_ref`) → reference RDKit index.
- A dock atom counts as the same atom iff OTMol maps it to the corresponding reference ligand atom. Written as the reference RDKit integer index (e.g. `L[7|8]`).
- **Hydrogens ignored**: RDKit ligand mols are loaded with `removeHs=True`, so Hs have no RDKit/OTMol counterpart; signatures keep heavy atoms only (H-bonds/salt bridges pair on the heavy donor/acceptor).
- **Ref coordinate map (key fix)**: BINANA relabels HETATM ligands to chain `X`/`UNK` while RDKit `mol_ref` is chain `B`, so `(chain,resID,resName,atomName)` tuple matching fails — the old ref signature silently fell back to raw PDB serials, which never matched the dock RDKit-index signatures, collapsing BINANA recall to ~0 for small-molecule ligands (and understating it for peptides). Ref now uses a coordinate map to RDKit indices just like dock, so both poses share one index space.

**Receptor part (the ONLY difference between levels; no OTMol)**:

| level | column | receptor label | granularity |
|---|---|---|---|
| eDockQ6_1 | default | chain : resID : resName : atomName | atom |
| eDockQ6_2 | rxatm | chain : resID : resName : atomName | atom (currently == 6_1) |
| eDockQ6_3 | rxidx | chain : resID : resName(upper) : atomName | atom, case-insensitive |
| eDockQ6_4 | rxres | chain : resID : resName (drops atomName) | residue (loosest) |

**Level-difference example (1B9J model 45, hydrophobic)**: ligand atom 15 ↔ GLU32, ref touches CB/CG, dock touches CD:
- atom-level 6_1/6_2/6_3: `L[15]R[A:32:GLU:CB]` ≠ `L[15]R[A:32:GLU:CD]` → no match.
- residue-level 6_4: `L[15]R[A:32:GLU]` on both → match.
That level's F1 rises 0.324 → 0.378.

### 2d. term_4 = (PoseBusters 8-check pass fraction)²
bond_lengths, bond_angles, aromatic_ring_flatness, double_bond_flatness, internal_steric_clash, energy_ratio≤100, minimum_distance_to_protein>0.75, volume_overlap_with_protein.

### 2e. TM = protein TMscore

> Example 1B9J model 0 (term_4=1.0, TM=1.0 are placeholders): pooled counts Σtp=6, Σref=6, Σdock=7 → B6 micro-F1=12/13=0.9231 (hbond F1=0.909, salt F1=1.0); H(hydrophobic)=0.4865 → B6·H=0.4491. All four levels: `1.0 × (0.4491 + 0.8856 + 0.9882)/3 × 1.0 = 0.7743`.
>
> Example T1124 model 0 (small-molecule ligand, placeholder term_4/TM): ligand RMSD 0.45 but interface RMSD 2.28 — the ligand is translated as a whole relative to the protein, so most native interactions are not reproduced. hbond 1/1/7/1 (recall=1 but 6 FPs), hydrophobic F1=0.1587, B6 micro-F1=0.1818 → B6·H=0.0289, eDockQ6 ≈ `1.0 × (0.0289 + 0.302 + 0.9972)/3 = 0.4427`.

---

## 3. Implementation notes (ver5, refined)

- Ligand OTMol mapping is computed **once** per (case, model) in `runner.py`; the same `dock_to_ref` feeds all downstream ligand work.
- **No MCS**: no `FindMCS / rdFMCS / GetSubstructMatch` in code, only "no MCS" clarifying comments.
- **No silent skipping**: failed case/models are never faked; they are recorded in `e-DockQ_results_failures.csv` / `posebusters_results_failures.csv` with reasons.
- `score_edockq_stepwise` produces only eDockQ1 and eDockQ6_1..4; eDockQ4/5/7–15 branches were removed. B6 uses the six-type micro pooling (`add_binana_six_pool_columns`), not a product.

## 4. Reproduction scripts
- `scripts/run_one_case_terms.py` — all intermediate terms + final scores for one case.
- `scripts/show_level_signatures.py [CASE MODEL]` — BINANA signatures and a summary table per level.
