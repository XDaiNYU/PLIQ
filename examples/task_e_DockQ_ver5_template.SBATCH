#!/bin/bash
#SBATCH --job-name=eDockQ_pliq
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --account=torch_pr_149_chemistry

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"

# Directory that contains run_pliq_stepwise_hpc.py AND lets `import pliq` work
# (i.e. the pliq_ver5_refined distribution root, or any dir on PYTHONPATH after
# `pip install .`). Edit to your cluster path.
PLIQ_RUN_ROOT=/scratch/xd638/vast_xd638/2025_eDockQ/Posebusters/pliq_ver5_refined
cd "${PLIQ_RUN_ROOT}"

singularity exec --overlay /scratch/xd638/singularity_common/common.ext3:ro \
  /share/apps/images/cuda11.8.86-cudnn8.7-devel-ubuntu22.04.2.sif \
  bash -lc "set -euo pipefail
cd \"${SLURM_SUBMIT_DIR}\"
source /ext3/env2.sh
conda activate pliq
cd \"${PLIQ_RUN_ROOT}\"

python examples/run_pliq_stepwise_hpc.py \
  --csv-file /scratch/xd638/vast_xd638/2025_eDockQ/Posebusters/processed_summary_v3.csv \
  --slice-start 0 --slice-stop 5 \
  --path-ref-lig /scratch/xd638/vast_xd638/2025_eDockQ/Posebusters_Yaowen_boltz/e_DockQ/Yaowen/posebusters_paper_data/posebusters_benchmark_set \
  --path-ref-pro /scratch/xd638/vast_xd638/2025_eDockQ/Posebusters_Yaowen_boltz/e_DockQ/Yaowen/posebusters_paper_data/posebusters_benchmark_set \
  --path-docked /scratch/xd638/vast_xd638/Yaowen_Boltz2_Xuhang \
  --saved-csv \"${SLURM_SUBMIT_DIR}/e-DockQ_stepwise_01.csv\" \
  --per-model-dir \"${SLURM_SUBMIT_DIR}/per_model\" \
  --n-models 20 \
  --binana-obabel-ph 7.4
"
