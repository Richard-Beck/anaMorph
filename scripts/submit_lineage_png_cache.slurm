#!/bin/bash
#SBATCH --job-name=lineage_png
#SBATCH --output=logs/lineage_png_%A_%a.out
#SBATCH --error=logs/lineage_png_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --qos=normal
#SBATCH --array=1-1000

set -euo pipefail

ALL_RESULTS_DIR=${1:-"data/lineage_area"}
MANIFEST_CSV=${2:-"manifests/image_manifest.csv"}
OUT_ROOT=${3:-"data/lineage_png"}
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
REPO_ROOT=$(pwd)

mkdir -p logs

readarray -t JOB_FIELDS < <(
python - <<'PY' "$ALL_RESULTS_DIR" "$MANIFEST_CSV" "$OUT_ROOT" "$TASK_ID"
import csv
import glob
import os
import sys

all_results_dir, manifest_csv, out_root, task_id = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

manifest = {}
with open(manifest_csv, newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        if row.get("artifact_type") == "raw_image" and row.get("filename") and row.get("filepath"):
            manifest[row["filename"]] = row["filepath"]

jobs = []
for image_csv in sorted(glob.glob(os.path.join(all_results_dir, "lineage_*", "image_summary.csv"))):
    lineage_dir = os.path.basename(os.path.dirname(image_csv))
    lineage_id = lineage_dir.removeprefix("lineage_")
    with open(image_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        seen = set()
        for row in reader:
            filename = row.get("filename")
            if not filename or filename in seen or filename not in manifest:
                continue
            seen.add(filename)
            raw_path = manifest[filename]
            png_name = os.path.splitext(os.path.basename(filename))[0] + ".png"
            out_path = os.path.join(out_root, lineage_dir, png_name)
            jobs.append((lineage_id, raw_path, out_path))

if task_id < 1 or task_id > len(jobs):
    print("__OUT_OF_RANGE__")
    sys.exit(0)

job = jobs[task_id - 1]
for field in job:
    print(field)
PY
)

if [[ ${JOB_FIELDS[0]} == "__OUT_OF_RANGE__" ]]; then
  echo "Array task ${TASK_ID} is out of range for the available lineage images. Exiting."
  exit 0
fi

LINEAGE_ID=${JOB_FIELDS[0]}
RAW_PATH=${JOB_FIELDS[1]}
OUT_PATH=${JOB_FIELDS[2]}

echo "=========================================================="
echo "Job started at $(date)"
echo "Task ID: ${TASK_ID}"
echo "Lineage ID: ${LINEAGE_ID}"
echo "Raw path: ${RAW_PATH}"
echo "Out path: ${OUT_PATH}"
echo "=========================================================="

python scripts/cache_lineage_png.py \
  --raw-path "${RAW_PATH}" \
  --out-path "${OUT_PATH}"

echo "=========================================================="
echo "Finished at $(date)"
echo "=========================================================="
