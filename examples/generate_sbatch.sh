#!/usr/bin/env bash
# =============================================================================
# PLIQ ver5 SBATCH generator (port of pliq_ver4/generate_casp15_task.sh).
#
# Slices each dataset listed in run_jobs.csv into ROWS_PER_SLICE-sized chunks and
# writes one SBATCH per slice (patched from task_e_DockQ_ver5_template.SBATCH),
# plus submit_all.sh and .csv_index.tsv.
#
# Layout on HPC:
#   <root>/
#     pliq_ver5_refined/                 <- this package (PLIQ_RUN_ROOT, `import pliq`)
#     <DATA_UPLOAD_DIR>/run_jobs.csv     <- datasets to score
#
# Usage:
#   cd .../pliq_ver5_refined/examples
#   ./generate_sbatch.sh
#   bash submit_all.sh
#
# Only one dataset:
#   JOB_NAME=alphafold3_casp15_raw ./generate_sbatch.sh
#
# Env vars:
#   HPC_ROOT          default /scratch/xd638/vast_xd638/2025_eDockQ/posebench_CASP15
#   PLIQ_RUN_ROOT     python run dir (dir where `import pliq` works)
#   LOCAL_POSEBENCH   local mirror prefix to rewrite into HPC paths
#   DATA_UPLOAD_DIR   data package folder name (contains run_jobs.csv)
#   ROWS_PER_SLICE    default 5
#   N_MODELS          default 20
#   CLEAN_OLD=1       remove old tasks_* before generating (default 1)
#   NO_TMSCORE=1      add --no-tmscore (default: TM-score on)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${HPC_WORK_DIR:-$SCRIPT_DIR}"

HPC_ROOT="${HPC_ROOT:-/scratch/xd638/vast_xd638/2025_eDockQ/posebench_CASP15}"
PLIQ_RUN_ROOT="${PLIQ_RUN_ROOT:-/scratch/xd638/vast_xd638/2025_eDockQ/Posebusters/pliq_ver5_refined}"
LOCAL_POSEBENCH="${LOCAL_POSEBENCH:-/Users/xuhangdai/Desktop/Cursor_working_folder/eDockQ/posebench}"
DATA_UPLOAD_DIR="${DATA_UPLOAD_DIR:-hpc_upload_pliq_casp15_20260608}"

RUN_JOBS_CSV="${RUN_JOBS_CSV:-${HPC_ROOT}/${DATA_UPLOAD_DIR}/run_jobs.csv}"
TEMPLATE_FILE="${TEMPLATE_FILE:-task_e_DockQ_ver5_template.SBATCH}"
ROWS_PER_SLICE="${ROWS_PER_SLICE:-5}"
N_MODELS="${N_MODELS:-20}"
JOB_NAME_FILTER="${JOB_NAME:-}"
CLEAN_OLD="${CLEAN_OLD:-1}"
NO_TMSCORE="${NO_TMSCORE:-0}"

if [[ ! -f "$RUN_JOBS_CSV" ]]; then
  alt="../${DATA_UPLOAD_DIR}/run_jobs.csv"
  if [[ -f "$alt" ]]; then
    RUN_JOBS_CSV="$(cd "$(dirname "$alt")" && pwd)/$(basename "$alt")"
  elif [[ -f "${LOCAL_POSEBENCH}/${DATA_UPLOAD_DIR}/run_jobs.csv" ]]; then
    RUN_JOBS_CSV="${LOCAL_POSEBENCH}/${DATA_UPLOAD_DIR}/run_jobs.csv"
  else
    echo "Cannot find run_jobs.csv" >&2
    echo "  tried: ${HPC_ROOT}/${DATA_UPLOAD_DIR}/run_jobs.csv" >&2
    echo "         ../${DATA_UPLOAD_DIR}/run_jobs.csv" >&2
    exit 1
  fi
fi

if [[ ! -f "$TEMPLATE_FILE" ]]; then
  echo "Cannot find template: $TEMPLATE_FILE" >&2
  exit 1
fi

echo "=============================================================================="
echo "PLIQ ver5 SBATCH generator"
echo "  work dir:        $(pwd)"
echo "  HPC root:        ${HPC_ROOT}"
echo "  PLIQ_RUN_ROOT:   ${PLIQ_RUN_ROOT}"
echo "  data package:    ${DATA_UPLOAD_DIR}"
echo "  slice:           ${ROWS_PER_SLICE}  n-models: ${N_MODELS}"
echo "  TM-score:        $([[ "${NO_TMSCORE}" == 1 ]] && echo skip || echo on)"
echo "  run_jobs:        ${RUN_JOBS_CSV}"
[[ -n "$JOB_NAME_FILTER" ]] && echo "  only:            ${JOB_NAME_FILTER}"
echo "=============================================================================="

export HPC_ROOT PLIQ_RUN_ROOT LOCAL_POSEBENCH DATA_UPLOAD_DIR
export RUN_JOBS_CSV TEMPLATE_FILE ROWS_PER_SLICE N_MODELS JOB_NAME_FILTER CLEAN_OLD NO_TMSCORE

python3 <<'PY'
import csv
import os
import re
import shutil
import sys
from pathlib import Path

HPC_ROOT = os.environ["HPC_ROOT"].rstrip("/")
PLIQ_RUN_ROOT = os.environ["PLIQ_RUN_ROOT"].rstrip("/")
LOCAL_POSEBENCH = os.environ["LOCAL_POSEBENCH"].rstrip("/")
DATA_UPLOAD_DIR = os.environ["DATA_UPLOAD_DIR"]
RUN_JOBS_CSV = os.environ["RUN_JOBS_CSV"]
TEMPLATE_FILE = os.environ["TEMPLATE_FILE"]
ROWS_PER_SLICE = int(os.environ["ROWS_PER_SLICE"])
N_MODELS = int(os.environ["N_MODELS"])
JOB_NAME_FILTER = os.environ.get("JOB_NAME_FILTER", "").strip()
CLEAN_OLD = os.environ.get("CLEAN_OLD", "1") == "1"
NO_TMSCORE = os.environ.get("NO_TMSCORE", "0") == "1"


def to_hpc(p: str) -> str:
    p = (p or "").strip()
    if not p:
        return p
    if DATA_UPLOAD_DIR in p:
        idx = p.find(DATA_UPLOAD_DIR)
        return f"{HPC_ROOT}/{p[idx:]}"
    if p.startswith(LOCAL_POSEBENCH):
        suffix = p[len(LOCAL_POSEBENCH):].lstrip("/")
        return f"{HPC_ROOT}/{suffix}"
    return p


def patch_template(tpl: str, job_name: str, row: dict, s: int, e: int, out_csv: str, per_model: str) -> str:
    r = tpl
    r = re.sub(r"(#SBATCH\s+--job-name=)[^\n]+", rf"\1{job_name}_{s:03d}_{e:03d}", r)
    r = re.sub(r"(PLIQ_RUN_ROOT=)[^\n]+", rf"\1{PLIQ_RUN_ROOT}", r)

    r = re.sub(r"\s*--no-binana\s*\\?\n?\s*", "\n  ", r)
    r = re.sub(r"\s*--no-posebusters\s*\\?\n?\s*", "\n  ", r)
    if not NO_TMSCORE:
        r = re.sub(r"\s*--no-tmscore\s*\\?\n?\s*", "\n  ", r)

    r = re.sub(r"--csv-file\s+\S+", f"--csv-file {row['csv_file']}", r)
    r = re.sub(r"--path-ref-lig\s+\S+", f"--path-ref-lig {row['path_ref_lig']}", r)
    r = re.sub(r"--path-ref-pro\s+\S+", f"--path-ref-pro {row['path_ref_pro']}", r)
    r = re.sub(r"--path-docked\s+\S+", f"--path-docked {row['path_docked']}", r)
    r = re.sub(r"--slice-start\s+\d+", f"--slice-start {s}", r)
    r = re.sub(r"--slice-stop\s+\d+", f"--slice-stop {e}", r)

    saved = f'--saved-csv "${{SLURM_SUBMIT_DIR}}/{out_csv}"'
    r = re.sub(r'--saved-csv\s+(?:"[^"]*"|\\"[^\\"]*\\")', saved, r)

    per = f'--per-model-dir "${{SLURM_SUBMIT_DIR}}/{per_model}"'
    r = re.sub(r'--per-model-dir\s+(?:"[^"]*"|\\"[^\\"]*\\")', per, r)

    if NO_TMSCORE and "--no-tmscore" not in r:
        r = r.replace(saved, f"--no-tmscore \\\n  {saved}")

    if re.search(r"--n-models\s+\d+", r):
        r = re.sub(r"--n-models\s+\d+", f"--n-models {N_MODELS}", r)
    else:
        r = re.sub(r"(--slice-stop\s+\d+)", rf"\1 \\\n  --n-models {N_MODELS}", r, count=1)
    return r


tpl = Path(TEMPLATE_FILE).read_text(encoding="utf-8")

if CLEAN_OLD:
    for d in sorted(Path(".").glob("tasks_*")):
        if d.is_dir():
            shutil.rmtree(d)

submit_lines = [
    "#!/bin/bash",
    "set -euo pipefail",
    'echo "submitting PLIQ ver5 SBATCH..."',
    "",
]
total = 0
index_rows = ["job_name\tcsv_no\tslice_start\tslice_stop\tcsv_file\n"]

print("\nDatasets:\n")
print(f"{'job_name':<40} {'n_cases':>8} {'tsv':>6} {'slices':>8}")
print("-" * 70)

with Path(RUN_JOBS_CSV).open(newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        job_name = row["job_name"].strip()
        if JOB_NAME_FILTER and job_name != JOB_NAME_FILTER:
            continue

        n_cases = int(row["n_cases"])
        csv_file = to_hpc(row["csv_file"].strip())
        path_ref_lig = to_hpc(row["path_ref_lig"].strip())
        path_ref_pro = to_hpc(row["path_ref_pro"].strip())
        path_docked = to_hpc(row["path_docked"].strip())

        tsv_rows = None
        if Path(csv_file).is_file():
            with Path(csv_file).open(encoding="utf-8") as tf:
                tsv_rows = sum(1 for _ in tf) - 1
            if tsv_rows != n_cases:
                print(f"  {job_name}: n_cases={n_cases} vs TSV {tsv_rows}, using TSV")
                n_cases = tsv_rows

        n_slices = (n_cases + ROWS_PER_SLICE - 1) // ROWS_PER_SLICE if n_cases else 0
        print(f"{job_name:<40} {int(row['n_cases']):>8} {str(tsv_rows or '?'):>6} {n_slices:>8}")

        out_dir = Path(f"tasks_{job_name}")
        out_dir.mkdir(parents=True, exist_ok=True)

        row_data = {
            "csv_file": csv_file,
            "path_ref_lig": path_ref_lig,
            "path_ref_pro": path_ref_pro,
            "path_docked": path_docked,
        }

        c = 1
        s = 0
        while s < n_cases:
            e = min(s + ROWS_PER_SLICE, n_cases)
            out_csv = f"tasks_{job_name}/e-DockQ_stepwise_{c:02d}.csv"
            per_model = f"tasks_{job_name}/per_model_{c:02d}"
            index_rows.append(f"{job_name}\t{c:02d}\t{s}\t{e}\t{out_csv}\n")

            sbatch_path = out_dir / f"task_{job_name}_slice_{s:03d}_{e:03d}.SBATCH"
            sbatch_path.write_text(
                patch_template(tpl, job_name, row_data, s, e, out_csv, per_model),
                encoding="utf-8",
            )

            submit_lines.append(f'echo "sbatch {sbatch_path}"')
            submit_lines.append(f"sbatch {sbatch_path}")
            total += 1
            c += 1
            s = e

Path("submit_all.sh").write_text("\n".join(submit_lines) + "\n", encoding="utf-8")
os.chmod("submit_all.sh", 0o755)
Path(".csv_index.tsv").write_text("".join(index_rows), encoding="utf-8")

print("\n" + "=" * 70)
print(f"Generated {total} SBATCH files.")
print("Submit:  bash submit_all.sh")
print("Index:   .csv_index.tsv")
print("=" * 70)

# Refuse to write local Mac paths into SBATCH (they won't exist on HPC).
bad = []
for sb in Path(".").glob("tasks_*/*.SBATCH"):
    txt = sb.read_text(encoding="utf-8")
    if "/Users/" in txt or "Desktop/Cursor" in txt:
        bad.append(str(sb))
if bad:
    print("\nERROR: these SBATCH still contain local paths (re-run on HPC):", file=sys.stderr)
    for b in bad[:5]:
        print(f"   {b}", file=sys.stderr)
    if len(bad) > 5:
        print(f"   ... {len(bad)} total", file=sys.stderr)
    raise SystemExit(1)
PY
