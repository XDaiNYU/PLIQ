# Contributing

Thank you for your interest in PLIQ.

## Development setup

```bash
git clone <your-fork-url>
cd PLIQ
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev]"     # or: pip install -e .
```

Optional external tools (not installed via pip):

- **Open Babel** (`obabel`) — BINANA PDB→PDBQT and PoseBusters MOL2 conversion
- **OTMol** (`import otmol`) — ligand alignment (no MCS fallback)
- **C++ compiler** — builds bundled TM-score on first use

## Running tests

```bash
pytest tests/ -q
```

Some tests expect benchmark structures under `data_posebusters/` (not shipped).
Mapping / unit tests should run without full data.

## Pull requests

1. Fork the repository and create a feature branch.
2. Keep changes focused; match existing code style and naming.
3. Add or update tests when behavior changes.
4. Update `CHANGELOG.md` for user-visible changes.
5. Open a PR with a clear description and test plan.

## Reporting issues

Please include:

- OS and Python version (`pliq info` output)
- Minimal command to reproduce
- Relevant log lines or `*_failures.csv` excerpts
