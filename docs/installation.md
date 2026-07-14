# Installation

## Requirements

- Python **≥ 3.8**
- C++ compiler (optional but recommended): builds TM-score from bundled source
- **Open Babel** (`obabel`) on `PATH` for BINANA and PoseBusters MOL2 paths
- **OTMol** (`import otmol`) for ligand alignment — required for e-DockQ metrics

## Install from source

```bash
cd PLIQ
pip install .
# editable dev install:
pip install -e ".[dev]"
```

Verify:

```bash
pliq --version
python -c "import pliq; print(pliq.__version__)"
pliq info
```

If `pliq` is not found, use `python -m pip install .` and `python -m pliq info`.

## Python dependencies (via pip)

| Package | Role |
|---------|------|
| pandas, numpy, tqdm | Tables and progress |
| prody, rdkit | Structure I/O |
| posebusters | Pose validation |
| openpyxl | Optional BINANA criteria override |

## TM-score binary

PLIQ does **not** ship a pre-built TM-score executable.

Resolution order:

1. `PLIQ_TM_EXE` environment variable
2. `$PLIQ_EDOCKQ_ROOT/TMscore`
3. Compile bundled `pliq/_vendor/tmscore/TMscore.cpp` into the user cache

Pre-build manually:

```bash
pliq compile-tmscore
```

Cache locations:

- macOS: `~/Library/Caches/pliq/TMscore`
- Linux: `~/.cache/pliq/TMscore`

Override cache root: `PLIQ_CACHE`.

## Environment variables

| Variable | Purpose |
|----------|---------|
| `PLIQ_EDOCKQ_ROOT` | Root containing `data_posebusters/` (batch `run-all`) |
| `PLIQ_TM_EXE` | Path to TM-score binary |
| `PLIQ_CACHE` | Override cache directory |

## Data not included

Large benchmark structures are **not** in this repository:

- `data_posebusters/` — reference + docked poses per case
- `otmol_package/` — may sit next to data root and is auto-detected

For HPC stepwise runs you only need your own `references/` and `docked/` trees
(see [HPC guide](hpc.md)).
