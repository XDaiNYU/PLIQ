# Protein-Ligand Interaction Quality (PLIQ) v6

Python package for end-to-end docking pose QA: ligand RMSD, interface RMSD, FNAT, **BINANA**, **PoseBusters**, **TM-score**, and a merged **PLIQ** score.

Ligand atom correspondence is **OTMol-only**; **reflection is always disabled**. BINANA uses **built-in default cutoffs**. Hydrogen handling uses **scheme A by default** (no `obabel -p`); see [BINANA hydrogen modes](#binana-hydrogen-modes) below.

**Docs:** [tutorial](docs/PLIQ_tutorial.ipynb) · [scoring](docs/scoring.md) · [HPC](docs/hpc.md) · [docs index](docs/README.md)

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

## One-command tutorial

Four PDBs in → one `pliq` score out:

```python
from pliq import run_pliq_from_pdbs

df = run_pliq_from_pdbs(
    "ref_ligand.pdb", "ref_protein.pdb",
    "dock_ligand.pdb", "dock_protein.pdb",
    "pliq_result.csv",
)
print(df["pliq"].iloc[0])
```

To reproduce the bundled example after cloning this repo:

```bash
python examples/7ZU2_DHT_model15/run_example.py
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

Upstream repo: **[https://github.com/weixiaoqimath/otmol](https://github.com/weixiaoqimath/otmol)**

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

**Check:** `pliq info` must show `[OK] OTMol`. Without OTMol, ligand-mapping / BINANA rows fail (recorded in `*_failures.csv`).

---



## TM-score — required for the PLIQ score (not a pip package)

**You do not need to download or point to TM-score manually in the normal case.**

PLIQ ships `TMscore.cpp` inside the package (`pliq/_vendor/tmscore/TMscore.cpp`). On first TM-score run it compiles a local binary (needs `g++` / `c++` in `PATH`).

### Default workflow (no path needed)

```bash
pip install .
pliq compile-tmscore          # optional pre-build; otherwise built on first run
pliq info                    # should show TM-score binary OK
```

Compiled binary location (automatic):


| OS      | Default path                                                |
| ------- | ----------------------------------------------------------- |
| macOS   | `~/Library/Caches/pliq/TMscore`                             |
| Linux   | `~/.cache/pliq/TMscore` (or `$XDG_CACHE_HOME/pliq/TMscore`) |
| Windows | `%LOCALAPPDATA%\pliq\cache\TMscore.exe`                     |


Override cache root: `export PLIQ_CACHE=/your/cache/dir`

### Do I need to set a TM-score path?


| Situation                                           | Action                                                                                 |
| --------------------------------------------------- | -------------------------------------------------------------------------------------- |
| **Typical user**                                    | **Nothing.** Run `pliq compile-tmscore` or let the first pipeline run compile it.      |
| **You already have a TM-score binary**              | Optional: `export PLIQ_TM_EXE=/path/to/TMscore`                                        |
| **You want to compile from your own** `TMscore.cpp` | Optional: `export PLIQ_TM_CPP=/path/to/TMscore.cpp` (with `PLIQ_TM_EXE`)               |
| **Binary in the data root**                         | If a `TMscore` executable sits next to `data_posebusters/`, PLIQ uses it automatically |


Legacy env aliases still work: `PLQ_TM_EXE`, `PLQ_TM_CPP`, `PLQ_CACHE`.

The data-root folder is for PDB trees (`data_posebusters/`), not for TM-score. TM-score resolution is independent unless you place a `TMscore` binary in that folder.

Upstream reference: **[https://zhanggroup.org/TM-score/](https://zhanggroup.org/TM-score/)** (use according to their license).

---



## Other dependencies


| Component                        | PyPI? | Install                                                             |
| -------------------------------- | ----- | ------------------------------------------------------------------- |
| **PLIQ**                         | —     | `pip install .` in `pliq_ver6/`                                     |
| **Open Babel**                   | No    | `conda install -c conda-forge openbabel` (needs `obabel` on `PATH`) |
| **BINANA**                       | —     | **bundled** in `pliq/_vendor/binana_python`                         |
| **PoseBusters, RDKit, ProDy, …** | Yes   | installed by `pip install .`                                        |




### External dependencies — quick reference


| Component      | PyPI?  | Where to get it                                                          | Install command                                                                              |
| -------------- | ------ | ------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------- |
| **PLIQ**       | —      | this repo (`pliq_ver6/`)                                                 | `pip install .`                                                                              |
| **OTMol**      | **No** | [github.com/weixiaoqimath/otmol](https://github.com/weixiaoqimath/otmol) | included in `pip install .`; or `pip install git+https://github.com/weixiaoqimath/otmol.git` |
| **TM-score**   | **No** | bundled `TMscore.cpp`                                                    | `pliq compile-tmscore` — **no manual path required**                                         |
| **Open Babel** | **No** | system / conda                                                           | `conda install -c conda-forge openbabel`                                                     |
| **BINANA**     | —      | **bundled** in `pliq/_vendor/binana_python`                              | nothing extra                                                                                |




## Package layout

```text
pliq_ver6/
  pyproject.toml          # pip metadata; declares otmol git dependency
  README.md
  pliq/                   # Python package
    __init__.py           # run_pliq_from_pdbs()
    deps.py               # dependency URLs + pliq info checks
    _vendor/
      binana_python/      # BINANA (bundled)
      tmscore/TMscore.cpp # TM-score source (bundled; build locally)
  examples/
  tests/
```



## End-to-end API (single pose) — `run_pliq_from_pdbs`

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
# df columns include ligand_rmsd, interface_rmsd, fnat_*, binana_*_tp_n,
# pliq (headline score), pliq_term_*
```

**Pipeline inside** `run_pliq_from_pdbs`**:**

1. OTMol ligand mapping + BINANA (`evaluate_pose_from_pdbs`)
2. PoseBusters (redock config)
3. TM-score (protein)
4. `pliq` **/** `pliq_term_`* **summary columns**



### PLIQ score columns (in `df` after scoring)

**Default headline score:** `pliq` (BINANA receptor level: chain + resID + resName + atomName).


| Column | Meaning                |
| ------ | ---------------------- |
| `pliq` | **Default PLIQ score** |


**Term decomposition for** `pliq`**:**

```
pliq = pliq_term_posebusters × (pliq_term_B6xH + pliq_term_kernel_i + pliq_term_kernel_l) / 3 × pliq_term_TM
```


| Column                  | Term         | Source                                          |
| ----------------------- | ------------ | ----------------------------------------------- |
| `pliq_term_fnat`        | BBSC FNAT F1 | backbone ∪ sidechain contact F1                 |
| `pliq_term_kernel_i`    | k_i          | `1 / (1 + (interface_rmsd/1.5)²)`               |
| `pliq_term_kernel_l`    | k_l          | `1 / (1 + (ligand_rmsd/8.5)²)`                  |
| `pliq_term_B6`          | B6           | Six-type BINANA micro-pooled F1 (default level) |
| `pliq_term_H`           | H            | Hydrophobic BINANA F1 (default level)           |
| `pliq_term_B6xH`        | B6×H         | `pliq_term_B6 × pliq_term_H`                    |
| `pliq_term_posebusters` | term_4       | `(PoseBusters 8-check pass fraction)²`          |
| `pliq_term_TM`          | TM           | TM-score (`tm_TMscore`)                         |


**Run-status flags** (boolean, end-to-end only):


| Column               | Module                            |
| -------------------- | --------------------------------- |
| `run_otmol_ok`       | OTMol ligand mapping + geometry   |
| `run_binana_ok`      | BINANA recall (no `binana_error`) |
| `run_posebusters_ok` | PoseBusters redock checks         |
| `run_tmscore_ok`     | TM-score protein alignment        |


---



## Worked example: 7ZU2_DHT model 15 (on GitHub)

The repo ships a **complete mini case** under `[examples/7ZU2_DHT_model15/](examples/7ZU2_DHT_model15/)`:

- **4 PDB files** (ref ligand/protein + dock ligand/protein)
- `run_example.py` — one command to reproduce
- `pliq_result_summary.csv` — 29 key columns with expected values
- `pliq_result_full.csv` — full ~740-column pipeline output
- **[Column guide](examples/7ZU2_DHT_model15/README.md)** — every summary column explained

```bash
pip install git+https://github.com/XDaiNYU/PLIQ.git
pliq compile-tmscore
cd examples/7ZU2_DHT_model15
python run_example.py
```

Expected headline score for this pose: `pliq` **≈ 0.928** (all four modules `run_*_ok = True`).

**Options:**


| Parameter         | Default | Meaning                                                 |
| ----------------- | ------- | ------------------------------------------------------- |
| `binana_strip_h`  | `False` | Scheme A: keep dock H. `True` = scheme B (`obabel -d`). |
| `run_posebusters` | `True`  | PoseBusters checks                                      |
| `run_tmscore`     | `True`  | TM-score                                                |
| `score_pliq`      | `True`  | Compute the PLIQ score and `pliq_term_*` columns        |


**Scheme B example:**

```python
df = run_pliq_from_pdbs(..., binana_strip_h=True)
```

**Metrics only (no PoseBusters / TM-score):**

```python
from pliq import evaluate_pose_from_pdbs

row = evaluate_pose_from_pdbs(ref_lig, ref_pro, dock_lig, dock_pro)
```

Defaults: OTMol alignment (no reflection), BINANA built-in cutoffs, PoseBusters redock, TM-score, and PLIQ scoring. Disable steps with `run_posebusters=False`, `run_tmscore=False`, or `score_pliq=False`.

## BINANA hydrogen modes

Neither scheme uses `obabel -p` (no added polar H).


| Mode            | `binana_strip_h` | Behavior                                                                          |
| --------------- | ---------------- | --------------------------------------------------------------------------------- |
| **A (default)** | `False`          | Ref PDBs as-is (crystal: no H). Dock keeps predictor protein/ligand H if present. |
| **B**           | `True`           | `obabel -d` strips H from **ref + dock** protein **and** ligand before BINANA.    |


CLI / batch: `--binana-strip-h` for scheme B.

Output column `binana_h_mode`: `as_is` or `strip_all`.

**Note:** Small-molecule dock ligands are typically heavy-atom only; dock **protein** H (e.g. AF3) is what matters for scheme A. Peptide ligands (NCAA) may carry H on both protein and ligand.

## What `pip install` delivers


| Component             | How it is delivered                                                                               |
| --------------------- | ------------------------------------------------------------------------------------------------- |
| **PLIQ Python code**  | Installed as the `pliq` package                                                                   |
| **PyPI dependencies** | `pandas`, `numpy`, `tqdm`, `prody`, `rdkit`, `posebusters`, `openpyxl`                            |
| **OTMol**             | Installed from [weixiaoqimath/otmol](https://github.com/weixiaoqimath/otmol) via `pyproject.toml` |
| **BINANA (Python)**   | **Bundled** under `pliq/_vendor/binana_python`                                                    |
| **TM-score**          | **Source** bundled; **binary built locally** (`pliq compile-tmscore`)                             |




## TM-score binary resolution (advanced)

Only needed if you use a **custom** binary instead of the default auto-build:

1. `PLIQ_TM_EXE` — if set and the file exists → PLIQ uses this binary (optional).
2. `TMscore` **in the data root** — if executable → PLIQ uses it (optional).
3. **Default** — compile bundled `pliq/_vendor/tmscore/TMscore.cpp` into the user cache (`pliq compile-tmscore` or on first run).

Optional env vars:

```bash
export PLIQ_TM_EXE=/path/to/TMscore      # optional: your pre-built binary
export PLIQ_TM_CPP=/path/to/TMscore.cpp   # optional: source for custom build
export PLIQ_CACHE=/path/to/cache          # optional: where to put compiled binary
```

Most users: **skip the above** and run `pliq compile-tmscore` once.

## Structures and system tools (not shipped)


| Need                                      | Why                                       |
| ----------------------------------------- | ----------------------------------------- |
| `data_posebusters/` or your own PDB trees | Input structures; large, not in the wheel |
| `obabel` on `PATH`                        | BINANA PDB→PDBQT and PoseBusters MOL2     |


```bash
e# data root = the folder that contains data_posebusters/
pliq run-all --output-dir ./pliq_runs/out
```

For single-pose evaluation you only need the four PDB files; no `data_posebusters/` layout is required.

## Failure handling (no silent skipping)

Dependencies are never faked. If a case/model cannot be computed, it is **skipped but explicitly recorded** in `*_failures.csv` with the error reason.

If `import otmol` is unavailable, affected ligand-mapping rows fail and are listed in the failures CSV rather than producing fabricated metrics.

## CLI

```text
pliq run-all [--cases CASES.tsv] [--output-dir DIR] [--n-models N] [--no-binana]
pliq compile-tmscore
pliq info
python -m pliq run-all ...
```

Outputs: results CSVs plus `end2end_<pdb_id>_merged.csv` per case.

## HPC: stepwise runner (same workflow as previous versions)

For Slurm, use the per-case **stepwise** runner — identical CLI to the old
`run_plq_stepwise_hpc_ver*.py`, only `plq`→`pliq`. It takes explicit ref/docked paths
plus an iloc slice (so each SBATCH does a chunk of cases) and writes **raw** metrics; no
`data_posebusters/` is required.

```bash
python examples/run_pliq_stepwise_hpc.py \
  --csv-file CASES.csv --slice-start 0 --slice-stop 5 \
  --path-ref-lig REF_DIR --path-ref-pro REF_DIR --path-docked DOCKED_DIR \
  --saved-csv pliq_stepwise_01.csv --per-model-dir per_model \
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

Score the raw stepwise CSV afterwards so it has `pliq` / `pliq_term_*` columns (see [docs/scoring.md](docs/scoring.md)).

Slurm helpers in `examples/`:

- SBATCH job template — edit `PLIQ_RUN_ROOT`, the singularity image, `conda activate <env>`, and paths.
- `generate_sbatch.sh` — slices every dataset in a `run_jobs.csv` into one SBATCH per chunk and writes `submit_all.sh`.

> Note: `--no-tmscore` drops the `tm_TMscore` column, which makes `pliq` **NaN**
> (`pliq` = `term_4 × (B6×H + k_i + k_l)/3 × TM`). Leave TM-score on for a complete score. On the
> compute node TM-score is compiled from the bundled `TMscore.cpp` on first use.



## Scores

The merged raw-metrics table is scored as:

```text
pliq = term_4 × (B6×H + k_i + k_l) / 3 × TM
```

See [docs/scoring.md](docs/scoring.md) for the per-term derivation.

## TM-score license

`TMscore.cpp` is from the Zhang group TM-score distribution; use according to their terms.