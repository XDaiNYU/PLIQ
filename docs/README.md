# Documentation

| Document | Description |
|----------|-------------|
| [Tutorial notebook](PLIQ_tutorial.ipynb) | One-pose walkthrough (OTMol → BINANA → PLIQ score) |
| [Installation](installation.md) | `pip install`, dependencies, TM-score build |
| [HPC & stepwise runner](hpc.md) | Slurm, `run_pliq_stepwise_hpc.py`, input layout |
| [Scoring & terms](scoring.md) | `pliq` formula, OTMol mapping overview |
| [Third-party licenses](third_party.md) | TM-score, BINANA |
| [Worked example](../examples/7ZU2_DHT_model15/README.md) | 7ZU2 / DHT / model 15 column guide |

## One-command tutorial

```python
from pliq import run_pliq_from_pdbs
df = run_pliq_from_pdbs("ref_lig.pdb", "ref_pro.pdb", "dock_lig.pdb", "dock_pro.pdb", "pliq_result.csv")
print(df["pliq"].iloc[0])
```

Bundled example: `python examples/7ZU2_DHT_model15/run_example.py`

## Quick links

```bash
pliq --version
pliq info
pliq run-all --help
python examples/run_pliq_stepwise_hpc.py --help
```

Score a stepwise CSV so it has `pliq` / `pliq_term_*` columns — see [scoring.md](scoring.md).
