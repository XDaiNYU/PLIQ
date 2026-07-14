# Protein-Ligand Interaction Quality (PLIQ) v6

Python package for end-to-end docking pose QA: **e-DockQ** (ligand RMSD, interface RMSD, FNAT, **BINANA**), **PoseBusters**, **TM-score**, merged CSVs.

Ligand atom correspondence is **OTMol-only**; **reflection is always disabled**. BINANA uses **built-in default cutoffs**. Hydrogen handling uses **scheme A by default** (no `obabel -p`); see [BINANA hydrogen modes](#binana-hydrogen-modes) below.

## Installation (Python package)

Requires **Python 3.9+**, **git** (for OTMol), and a **C++ compiler** (`g++` or `c++`, for TM-score on first use).

```bash
cd pliq_ver6
pip install .
```

This installs the `pliq` package, PyPI dependencies, and **OTMol from GitHub** (declared in `pyproject.toml`).

Verify everything:

```bash
pliq --version
pliq info    # reports OTMol, obabel, TM-score binary status
```

---

## OTMol — required (ligand RMSD + BINANA mapping)

**OTMol is not on PyPI.** It must be importable as `import otmol`.

### Recommended: install with PLIQ (automatic)

`pip install .` already pulls OTMol from GitHub:

```bash
cd pliq_ver6
pip install .
python -c "import otmol; print('otmol OK')"
```

Upstream repo: **https://github.com/weixiaoqimath/otmol**

### Standalone install (if `import otmol` fails)

```bash
pip install git+https://github.com/weixiaoqimath/otmol.git
python -c "import otmol; print('otmol OK')"
```

### Alternative: conda env from upstream

```bash
git clone https://github.com/weixiaoqimath/otmol.git
cd otmol
conda env create -f environment.yml
conda activate otmol    # env name per upstream environment.yml
pip install .
```

**Check:** `pliq info` must show `[OK] OTMol`. Without OTMol, e-DockQ / BINANA rows fail (recorded in `*_failures.csv`).

---

## TM-score — required for eDockQ6 (not a pip package)

**You do not need to download or point to TM-score manually in the normal case.**

PLIQ ships `TMscore.cpp` inside the package (`pliq/_vendor/tmscore/TMscore.cpp`). On first TM-score run it compiles a local binary (needs `g++` / `c++` in `PATH`).

### Default workflow (no path needed)

```bash
pip install .
pliq compile-tmscore          # optional pre-build; otherwise built on first run
pliq info                    # should show TM-score binary OK
```

Compiled binary location (automatic):

| OS | Default path |
|----|----------------|
| macOS | `~/Library/Caches/pliq/TMscore` |
| Linux | `~/.cache/pliq/TMscore` (or `$XDG_CACHE_HOME/pliq/TMscore`) |
| Windows | `%LOCALAPPDATA%\pliq\cache\TMscore.exe` |

Override cache root: `export PLIQ_CACHE=/your/cache/dir`

### Do I need to set a TM-score path?

| Situation | Action |
|-----------|--------|
| **Typical user** | **Nothing.** Run `pliq compile-tmscore` or let the first pipeline run compile it. |
| **You already have a TM-score binary** | Optional: `export PLIQ_TM_EXE=/path/to/TMscore` |
| **You want to compile from your own `TMscore.cpp`** | Optional: `export PLIQ_TM_CPP=/path/to/TMscore.cpp` (with `PLIQ_TM_EXE`) |
| **Binary in data repo root** | If `$PLIQ_EDOCKQ_ROOT/TMscore` exists and is executable, PLIQ uses it automatically |

Legacy env aliases still work: `PLQ_TM_EXE`, `PLQ_TM_CPP`, `PLQ_CACHE`.

**`PLIQ_EDOCKQ_ROOT` is for PDB data (`data_posebusters/`), not for TM-score.** TM-score resolution is independent unless you place a `TMscore` binary at the repo root.

Upstream reference: **https://zhanggroup.org/TM-score/** (use according to their license).

---

## Other dependencies

| Component | PyPI? | Install |
|-----------|-------|---------|
| **PLIQ** | — | `pip install .` in `pliq_ver6/` |
| **Open Babel** | No | `conda install -c conda-forge openbabel` (needs `obabel` on `PATH`) |
| **BINANA** | — | **bundled** in `pliq/_vendor/binana_python` |
| **PoseBusters, RDKit, ProDy, …** | Yes | installed by `pip install .` |

### External dependencies — quick reference

| Component | PyPI? | Where to get it | Install command |
|-----------|-------|-----------------|-----------------|
| **PLIQ** | — | this repo (`pliq_ver6/`) | `pip install .` |
| **OTMol** | **No** | [github.com/weixiaoqimath/otmol](https://github.com/weixiaoqimath/otmol) | included in `pip install .`; or `pip install git+https://github.com/weixiaoqimath/otmol.git` |
| **TM-score** | **No** | bundled `TMscore.cpp` | `pliq compile-tmscore` — **no manual path required** |
| **Open Babel** | **No** | system / conda | `conda install -c conda-forge openbabel` |
| **BINANA** | — | **bundled** in `pliq/_vendor/binana_python` | nothing extra |

## Package layout

```text
pliq_ver6/
  pyproject.toml          # pip metadata; declares otmol git dependency
  README.md
  pliq_edockq_scoring.py
  pliq/                   # Python package
    __init__.py           # run_pliq_from_pdbs()
    deps.py               # dependency URLs + pliq info checks
    _vendor/
      binana_python/      # BINANA (bundled)
      tmscore/TMscore.cpp # TM-score source (bundled; build locally)
  examples/
  tests/
```

## End-to-end API (single pose) — **`run_pliq_from_pdbs`**

**Four PDB files in → one scored CSV out.** This is the ultimate single-pose entry point.

```python
from pliq import run_pliq_from_pdbs

df = run_pliq_from_pdbs(
    "ref_ligand.pdb",      # reference ligand (chain B)
    "ref_protein.pdb",     # reference protein (chain A)
    "dock_ligand.pdb",     # docked ligand (chain B)
    "dock_protein.pdb",    # docked protein (chain A)
    "pliq_result.csv",
    pdb_id="1ABC",
    af_model_id=0,
)
# df columns include ligand_rmsd, interface_rmsd, fnat_*, binana_*_tp_n, eDockQ1, eDockQ6_1..4
```

**Pipeline inside `run_pliq_from_pdbs`:**

1. OTMol e-DockQ + BINANA (`evaluate_pose_from_pdbs`)
2. PoseBusters (redock config)
3. TM-score (protein)
4. `eDockQ1` / `eDockQ6_1..4` scoring (`pliq_edockq_scoring`)

**Options:**

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `binana_strip_h` | `False` | Scheme A: keep dock H. `True` = scheme B (`obabel -d`). |
| `run_posebusters` | `True` | PoseBusters checks |
| `run_tmscore` | `True` | TM-score |
| `score_pliq` | `True` | Compute eDockQ1 / eDockQ6 |

**Scheme B example:**

```python
df = run_pliq_from_pdbs(..., binana_strip_h=True)
```

**Metrics only (no PoseBusters / TM-score):**

```python
from pliq import evaluate_pose_from_pdbs

row = evaluate_pose_from_pdbs(ref_lig, ref_pro, dock_lig, dock_pro)
```

Defaults: OTMol alignment (no reflection), BINANA built-in cutoffs, PoseBusters redock, TM-score, and `eDockQ1` / `eDockQ6_1..4` scoring. Disable steps with `run_posebusters=False`, `run_tmscore=False`, or `score_pliq=False`.

## BINANA hydrogen modes

Neither scheme uses `obabel -p` (no added polar H).

| Mode | `binana_strip_h` | Behavior |
|------|------------------|----------|
| **A (default)** | `False` | Ref PDBs as-is (crystal: no H). Dock keeps predictor protein/ligand H if present. |
| **B** | `True` | `obabel -d` strips H from **ref + dock** protein **and** ligand before BINANA. |

CLI / batch: `--binana-strip-h` for scheme B.

Output column `binana_h_mode`: `as_is` or `strip_all`.

**Note:** Small-molecule dock ligands are typically heavy-atom only; dock **protein** H (e.g. AF3) is what matters for scheme A. Peptide ligands (NCAA) may carry H on both protein and ligand.

## What `pip install` delivers

| Component | How it is delivered |
|-----------|---------------------|
| **PLIQ Python code** | Installed as the `pliq` package |
| **PyPI dependencies** | `pandas`, `numpy`, `tqdm`, `prody`, `rdkit`, `posebusters`, `openpyxl` |
| **OTMol** | Installed from [weixiaoqimath/otmol](https://github.com/weixiaoqimath/otmol) via `pyproject.toml` |
| **BINANA (Python)** | **Bundled** under `pliq/_vendor/binana_python` |
| **TM-score** | **Source** bundled; **binary built locally** (`pliq compile-tmscore`) |

## TM-score binary resolution (advanced)

Only needed if you use a **custom** binary instead of the default auto-build:

1. **`PLIQ_TM_EXE`** — if set and the file exists → PLIQ uses this binary (optional).
2. **`$PLIQ_EDOCKQ_ROOT/TMscore`** — if executable → PLIQ uses it (optional).
3. **Default** — compile bundled `pliq/_vendor/tmscore/TMscore.cpp` into the user cache (`pliq compile-tmscore` or on first run).

Optional env vars:

```bash
export PLIQ_TM_EXE=/path/to/TMscore      # optional: your pre-built binary
export PLIQ_TM_CPP=/path/to/TMscore.cpp   # optional: source for custom build
export PLIQ_CACHE=/path/to/cache          # optional: where to put compiled binary
```

Most users: **skip the above** and run `pliq compile-tmscore` once.

## Structures and system tools (not shipped)

| Need | Why |
|------|-----|
| **`data_posebusters/`** or your own PDB trees | Input structures; large, not in the wheel |
| **`obabel`** on `PATH` | BINANA PDB→PDBQT and PoseBusters MOL2 |

```bash
export PLIQ_EDOCKQ_ROOT=/path/to/eDockQ    # directory containing data_posebusters/
pliq run-all --output-dir ./pliq_runs/out
```

For single-pose evaluation you only need the four PDB files; no `data_posebusters/` layout is required.

## Failure handling (no silent skipping)

Dependencies are never faked. If a case/model cannot be computed, it is **skipped but explicitly recorded**:
- `e-DockQ_results_failures.csv` — per case/model with the error reason.
- `posebusters_results_failures.csv` — per case/model with the error reason.

If `import otmol` is unavailable, affected e-DockQ rows fail and are listed in the failures CSV rather than producing fabricated metrics.

## CLI

```text
pliq run-all [--cases CASES.tsv] [--edockq-root DIR] [--output-dir DIR] [--n-models N] [--no-binana]
pliq compile-tmscore
pliq info
python -m pliq run-all ...
```

Outputs: `e-DockQ_results.csv`, `posebusters_results.csv`, `TMscore_results.csv`, `merged_results.csv`, and `end2end_<pdb_id>_merged.csv` per case.

## HPC: stepwise runner (same workflow as previous versions)

For Slurm, use the per-case **stepwise** runner — identical CLI to the old
`run_plq_stepwise_hpc_ver*.py`, only `plq`→`pliq`. It takes explicit ref/docked paths
plus an iloc slice (so each SBATCH does a chunk of cases) and writes **raw** metrics; no
`data_posebusters/` is required.

```bash
python examples/run_pliq_stepwise_hpc.py \
  --csv-file CASES.csv --slice-start 0 --slice-stop 5 \
  --path-ref-lig REF_DIR --path-ref-pro REF_DIR --path-docked DOCKED_DIR \
  --saved-csv e-DockQ_stepwise_01.csv --per-model-dir per_model \
  --n-models 20
```

BINANA: scheme A by default (no `obabel -p`). Scheme B: add `--binana-strip-h`.

Per-case file layout it expects:

```text
{path_ref_lig}/{pdb_id}/{pdb_id}_ligand_chainB.pdb
{path_ref_pro}/{pdb_id}/{pdb_id}_protein_renum.pdb
{path_docked}/posebusters_{pdb_id}/posebusters_{pdb_id}_model_{k}_aligned_ligand.pdb
{path_docked}/posebusters_{pdb_id}/posebusters_{pdb_id}_model_{k}_aligned_protein.pdb
```

Turn the raw stepwise CSV into eDockQ scores afterwards:

```bash
python -c "import pandas as pd, pliq_edockq_scoring as S; \
  df=pd.read_csv('e-DockQ_stepwise_01.csv'); \
  S.score_edockq_stepwise(df).to_csv('e-DockQ_scored_01.csv', index=False)"
```

Slurm helpers in `examples/`:
- `task_e_DockQ_ver5_template.SBATCH` — edit `PLIQ_RUN_ROOT`, the singularity image, `conda activate <env>`, and paths.
- `generate_sbatch.sh` — slices every dataset in a `run_jobs.csv` into one SBATCH per chunk and writes `submit_all.sh` (port of the ver4 generator; `PLQ_*`→`PLIQ_*`).

> Note: `--no-tmscore` drops the `tm_TMscore` column, which makes **eDockQ6_1..4 NaN**
> (eDockQ6 = `term_4 × (B6×H + k_i + k_l)/3 × TM`). Leave TM-score on for eDockQ6. On the
> compute node TM-score is compiled from the bundled `TMscore.cpp` on first use.

## Scores

`pliq_edockq_scoring.score_edockq_stepwise(df)` turns the merged raw-metrics table into:
- **eDockQ1** = `(fnat_BBSC + k_i + k_l) / 3`
- **eDockQ6_1..4** = `term_4 × (B6×H + k_i + k_l)/3 × TM` for the four BINANA receptor-granularity levels.

See `reports/pliq_ver5_report_en.md` / `_zh.md` for the full per-term derivation and mapping details.

## TM-score license

`TMscore.cpp` is from the Zhang group TM-score distribution; use according to their terms.
