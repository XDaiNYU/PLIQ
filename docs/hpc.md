# HPC and stepwise evaluation

Use the **stepwise runner** when structures live on the cluster and you want
one SBATCH job per case slice — no bundled `data_posebusters/` required.

## Input layout

For each case `PDB_ID` (e.g. `7U0U_FK5`):

```text
{path_ref_lig}/{PDB_ID}/{PDB_ID}_ligand_chainB.pdb
{path_ref_pro}/{PDB_ID}/{PDB_ID}_protein_renum.pdb
{path_docked}/posebusters_{PDB_ID}/posebusters_{PDB_ID}_model_{k}_aligned_ligand.pdb
{path_docked}/posebusters_{PDB_ID}/posebusters_{PDB_ID}_model_{k}_aligned_protein.pdb
```

**Protein:** reference and dock should share **matched residue numbering**
(chain A, residues `1..M`, same `(resseq, resname)` pairs) for reliable
interface RMSD, FNAT, and TM-score.

**Ligand:** chain B; atom correspondence is **OTMol-only** (not resseq-based).

## Stepwise command

```bash
python examples/run_pliq_stepwise_hpc.py \
  --csv-file cases.tsv \
  --csv-sep $'\t' \
  --subfolder-col pdb_id \
  --slice-start 0 --slice-stop 1 \
  --path-ref-lig /path/references \
  --path-ref-pro /path/references \
  --path-docked /path/docked \
  --saved-csv pliq_stepwise_001.csv \
  --per-model-dir per_model_001 \
  --n-models 20 \
  --binana-obabel-ph 7.4
```

`cases.tsv` format:

```text
pdb_id	n_models
7U0U_FK5	20
```

## Slurm helpers

In `examples/`:

| File | Purpose |
|------|---------|
| SBATCH job template | Edit account, singularity, conda env |
| `generate_sbatch.sh` | Slice `run_jobs.csv` into per-chunk SBATCH + `submit_all.sh` |

On HPC:

```bash
export PLIQ_RUN_ROOT=/path/to/PLIQ   # this repo after pip install .
./generate_sbatch.sh
bash submit_all.sh
```

## Score the raw CSV

The stepwise runner writes **raw metrics**. Add `pliq` / `pliq_term_*` afterwards
(see [scoring.md](scoring.md)).

## TM-score

Leave TM-score **enabled** (default) so `pliq` is defined:

```text
pliq = term_4 × (B6×H + k_i + k_l)/3 × TM
```

Use `--no-tmscore` only if the binary cannot be built; the `pliq` column will be NaN.

## Outputs

Per job:

- `pliq_stepwise_*.csv` — merged raw metrics
- `per_model_*/*.csv` — one row per model
- Failures are recorded in-row (`*_error` columns) and printed at end
