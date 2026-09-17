# Scoring overview

PLIQ produces a headline **`pliq`** score (plus `pliq_term_*` columns) from a
single raw-metrics table.

## Pipeline

```text
Structures  →  run-all / stepwise runner  →  raw CSV
                                              ↓
                                    PLIQ scoring
                                              ↓
                                    pliq, pliq_term_*
```

For a single pose, `run_pliq_from_pdbs()` already writes the scored columns.

## PLIQ score

```text
pliq = term_4 × (B6×H + k_i + k_l) / 3 × TM
```

| Term | Column | Meaning |
|------|--------|---------|
| `fnat_BBSC` | `pliq_term_fnat` | Fraction native contacts (backbone + sidechain F1) |
| `k_i` | `pliq_term_kernel_i` | Interface RMSD kernel, `1/(1+(iRMSD/1.5)²)` |
| `k_l` | `pliq_term_kernel_l` | Ligand RMSD kernel, `1/(1+(ligRMSD/8.5)²)` |
| `B6` | `pliq_term_B6` | Micro-pooled F1 over six BINANA types (hbond, salt, pipi, tstack, cationpi, halogen) |
| `H` | `pliq_term_H` | Hydrophobic interaction F1 |
| `B6×H` | `pliq_term_B6xH` | `pliq_term_B6 × pliq_term_H` |
| `term_4` | `pliq_term_posebusters` | `(PoseBusters 8-check pass fraction)²` |
| `TM` | `pliq_term_TM` | Protein TM-score (ref vs dock) |

BINANA receptor granularity for the default `pliq` score is
chain + resID + resName + atomName.

## Ligand mapping (OTMol)

One `dock_to_ref` map per (case, model):

- Used by ligand RMSD, interface RMSD, FNAT, and all BINANA ligand signatures
- **No MCS** anywhere

Protein atoms are matched by `(resnum, resname, atom_name)` after inputs are
aligned to common numbering.

## Further reading

- [Tutorial notebook](PLIQ_tutorial.ipynb)
- [Worked example](../examples/7ZU2_DHT_model15/README.md)
- `scripts/run_one_case_terms.py` — print all terms for one case/model
