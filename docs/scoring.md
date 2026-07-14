# Scoring overview

PLIQ ver5 (refined) produces **eDockQ1** and **eDockQ6_1 .. eDockQ6_4** from a
single raw-metrics table.

## Pipeline

```text
Structures  →  run-all / stepwise runner  →  raw CSV
                                              ↓
                              pliq_edockq_scoring.score_edockq_stepwise()
                                              ↓
                                    eDockQ1, eDockQ6_1..4
```

## eDockQ1

```text
eDockQ1 = (fnat_BBSC + k_i + k_l) / 3
```

| Term | Meaning |
|------|---------|
| `fnat_BBSC` | Fraction native contacts (backbone + sidechain F1) |
| `k_i` | Interface RMSD kernel, `1/(1+(iRMSD/1.5)²)` |
| `k_l` | Ligand RMSD kernel, `1/(1+(ligRMSD/8.5)²)` |

## eDockQ6_k (k = 1..4)

```text
eDockQ6_k = term_4 × (B6×H + k_i + k_l)/3 × TM
```

| Term | Meaning |
|------|---------|
| `term_4` | PoseBusters-style validity gate |
| `B6` | Micro-pooled F1 over six BINANA types (hbond, salt, pipi, tstack, cationpi, halogen) |
| `H` | Hydrophobic interaction F1 |
| `TM` | Protein TM-score (ref vs dock) |

The four **k** levels differ only in **BINANA receptor granularity**
(atom / residue / index / rxatm modes).

## Ligand mapping (OTMol)

One `dock_to_ref` map per (case, model):

- Used by ligand RMSD, interface RMSD, FNAT, and all BINANA ligand signatures
- **No MCS** anywhere

Protein atoms are matched by `(resnum, resname, atom_name)` after inputs are
aligned to common numbering.

## Further reading

- [Full English report](../reports/pliq_ver5_report_en.md)
- [Full Chinese report](../reports/pliq_ver5_report_zh.md)
- `scripts/run_one_case_terms.py` — print all terms for one case/model
