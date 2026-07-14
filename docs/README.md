# Documentation

| Document | Description |
|----------|-------------|
| [Installation](installation.md) | `pip install`, dependencies, TM-score build |
| [HPC & stepwise runner](hpc.md) | Slurm, `run_pliq_stepwise_hpc.py`, input layout |
| [Scoring & terms](scoring.md) | eDockQ1, eDockQ6_1..4, OTMol mapping overview |
| [Third-party licenses](third_party.md) | TM-score, BINANA |
| [Full term report (EN)](../reports/pliq_ver5_report_en.md) | Detailed per-term derivation |
| [Full term report (ZH)](../reports/pliq_ver5_report_zh.md) | 中文详细说明 |

## Quick links

```bash
pliq --version
pliq info
pliq run-all --help
python examples/run_pliq_stepwise_hpc.py --help
```

Post-process stepwise CSV:

```bash
python -c "import pandas as pd, pliq_edockq_scoring as S; \
  df=pd.read_csv('e-DockQ_stepwise.csv'); \
  S.score_edockq_stepwise(df).to_csv('e-DockQ_scored.csv', index=False)"
```
