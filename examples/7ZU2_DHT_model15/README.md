# Worked example: 7ZU2_DHT model 15

Bundled reference + docked PDBs and pre-computed PLIQ output for a single PoseBuster pose.

## Files

```text
examples/7ZU2_DHT_model15/
  ref/
    7ZU2_DHT_ligand_chainB.pdb          # reference ligand (chain B)
    7ZU2_DHT_protein_renum.pdb          # reference protein (chain A)
  dock/
    posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb
    posebusters_7ZU2_DHT_model_15_aligned_protein.pdb
  run_example.py                        # reproduce the run
  pliq_result_summary.csv               # key columns (29 cols)
  pliq_result_full.csv                  # full pipeline output (~740 cols)
```

## How to run

```bash
pip install git+https://github.com/XDaiNYU/PLIQ.git
pliq compile-tmscore    # once, needs g++ / c++
pliq info             # OTMol + obabel + TM-score OK

cd examples/7ZU2_DHT_model15
python run_example.py
```

Or from Python anywhere:

```python
from pathlib import Path
from pliq import run_pliq_from_pdbs

root = Path("examples/7ZU2_DHT_model15")  # path after cloning PLIQ

df = run_pliq_from_pdbs(
    root / "ref/7ZU2_DHT_ligand_chainB.pdb",
    root / "ref/7ZU2_DHT_protein_renum.pdb",
    root / "dock/posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb",
    root / "dock/posebusters_7ZU2_DHT_model_15_aligned_protein.pdb",
    root / "pliq_result_full.csv",
    pdb_id="7ZU2_DHT",
    af_model_id=15,
)
print(df[["pliq", "run_otmol_ok", "run_binana_ok"]].iloc[0])
```

`run_pliq_from_pdbs` returns a **one-row `pandas.DataFrame`** and writes the same row to CSV.

## Expected result (this case)

| Column | Value |
|--------|------:|
| `run_otmol_ok` | True |
| `run_binana_ok` | True |
| `run_posebusters_ok` | True |
| `run_tmscore_ok` | True |
| `ligand_rmsd` | 0.127 |
| `interface_rmsd` | 0.193 |
| `fnat_all` | 0.935 |
| `tm_TMscore` | 0.992 |
| **`pliq`** | **0.928** |

## Column guide (`pliq_result_summary.csv`)

### Run status

| Column | Module | Meaning |
|--------|--------|---------|
| `run_otmol_ok` | OTMol | Ligand mapping + geometry succeeded |
| `run_binana_ok` | BINANA | Interaction recall ran (`binana_error` empty) |
| `run_posebusters_ok` | PoseBusters | Redock checks succeeded |
| `run_tmscore_ok` | TM-score | Protein TM-score computed |

### Raw metrics (selected)

| Column | Meaning |
|--------|---------|
| `pdb_id` / `af_model_id` | Case id and model index |
| `ligand_rmsd` | OTMol ligand RMSD (Å) |
| `interface_rmsd` | Interface CA RMSD (Å) |
| `fnat_all` | FNAT over all atom contacts |
| `otmol_map_pairs` | OTMol dock→ref atom pairs |
| `binana_h_mode` | `as_is` (scheme A) or `strip_all` (scheme B) |
| `binana_hbond_tp_n` | BINANA hydrogen-bond true positives |
| `binana_hydrophobic_tp_n` | BINANA hydrophobic true positives |
| `tm_TMscore` | TM-score |
| `pb_valid_fraction` | PoseBusters checks passed / 8 |

### PLIQ scores

| Column | Meaning |
|--------|---------|
| **`pliq`** | **Default headline score** |

### PLIQ term decomposition

```
pliq = pliq_term_posebusters × (pliq_term_B6xH + pliq_term_kernel_i + pliq_term_kernel_l) / 3 × pliq_term_TM
```

| Column | Term | Formula / source |
|--------|------|------------------|
| `pliq_term_fnat` | BBSC FNAT | Backbone ∪ sidechain contact F1 |
| `pliq_term_kernel_i` | k_i | `1 / (1 + (interface_rmsd/1.5)²)` |
| `pliq_term_kernel_l` | k_l | `1 / (1 + (ligand_rmsd/8.5)²)` |
| `pliq_term_B6` | B6 | Six-type BINANA micro-pooled F1 (default level) |
| `pliq_term_H` | H | Hydrophobic BINANA F1 (default level) |
| `pliq_term_B6xH` | B6×H | `pliq_term_B6 × pliq_term_H` |
| `pliq_term_posebusters` | term_4 | `(pb_valid_fraction)²` |
| `pliq_term_TM` | TM | `tm_TMscore` |

`pliq_result_full.csv` adds all intermediate BINANA detail columns, PoseBusters per-check columns (`pb_*`), and per-level term columns (`term_binana_*_{default,rxatm,rxidx,rxres}`).
