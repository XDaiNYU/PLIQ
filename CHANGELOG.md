# Changelog

All notable changes to **PLIQ** are documented here.

Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.6.1] — 2026-06-23

### Added

- **`run_pliq_from_pdbs()`** — ultimate end-to-end API: four PDB files → scored CSV (e-DockQ, BINANA, PoseBusters, TM-score, eDockQ1/eDockQ6).
- **`evaluate_pose_from_pdbs()`** — single-pose metrics without PoseBusters/TM-score.
- **BINANA hydrogen scheme B** — `binana_strip_h=True` / `--binana-strip-h` (`obabel -d` on ref+dock protein+ligand).
- **`pliq/deps.py`** and enhanced **`pliq info`** for OTMol / obabel / TM-score checks.
- **`tests/test_binana_strip_h.py`** for scheme A vs B.

### Changed

- **OTMol** installed from [weixiaoqimath/otmol](https://github.com/weixiaoqimath/otmol) via `pyproject.toml`.
- **Reflection always disabled** in ligand alignment presets.
- **BINANA default hydrogen (scheme A)**: no `obabel -p`; ref as-is (no crystal H), dock keeps predictor H if present.
- **BINANA built-in default cutoffs** unless optional criteria file is passed.
- README: OTMol / TM-score install guide; TM-score path not required by default.

## [0.5.0] — 2026-06-10

### Added

- **BINANA per-interaction TP / FN / FP details** with JSON detail columns and signature lists.
- **Coordinate-based BINANA → RDKit maps** for both reference and dock ligands.
- **`pliq_edockq_scoring.py`**: post-processing for eDockQ1 and eDockQ6_1..4 only.
- **HPC stepwise runner** `examples/run_pliq_stepwise_hpc.py` with Slurm helpers.
- **Failure CSVs** (`*_failures.csv`) instead of silent skips when dependencies fail.
- **`add_binana_six_pool_columns`**: B6 as one micro-pooled F1 over six interaction types.

### Changed

- **Single OTMol ligand map** per (case, model), reused by RMSD, FNAT, and all BINANA levels.
- **Four BINANA levels** differ only by receptor granularity; ligand side always uses OTMol indices.
- **Canonical chains** on read: ligand → `B`, protein → `A`.
- **Package rename** `plq` → `pliq`; removed legacy duplicate package and MCS paths.

### Fixed

- BINANA ref-ligand signature mismatch (chain `X` vs `B`) that collapsed recall to ~0.
- Ligand hydrogen atoms excluded from BINANA signatures (no OTMol counterpart).
- B6 formula: micro-pooled F1, not product of six F1s.

### Removed

- eDockQ4/5/7–15 and pose6 micro-F1 branches from scoring.
- Unused helpers and legacy `run_plq_stepwise*.py` examples.
