#!/bin/bash
#SBATCH --job-name=lineage_area
#SBATCH --output=logs/lineage_area_%A_%a.out
#SBATCH --error=logs/lineage_area_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --qos=normal
#SBATCH --array=1-88

set -euo pipefail

LINEAGE_RDS=${1:-"core_data/lineages.Rds"}
OBJECT_SIZES_RDS=${2:-"data/all_images/object_sizes.rds"}
OUT_DIR=${3:-"data/lineage_area"}
CONTAINER_URI=${4:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}

REPO_ROOT=$(pwd)
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
BINDS="-B ${REPO_ROOT},/home/$USER,/share,/etc/passwd,/etc/group"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export TZ="US/Eastern"
export APPTAINER_NOHTTPS=1

mkdir -p logs
mkdir -p "${OUT_DIR}"

echo "=========================================================="
echo "Job started at $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task: ${TASK_ID}"
echo "Repo root: ${REPO_ROOT}"
echo "Lineage RDS: ${LINEAGE_RDS}"
echo "Object sizes RDS: ${OBJECT_SIZES_RDS}"
echo "Output dir: ${OUT_DIR}"
echo "Container: ${CONTAINER_URI}"
echo "=========================================================="

LINEAGE_ID=$(apptainer exec ${BINDS} "${CONTAINER_URI}" \
  Rscript -e "x <- readRDS('${LINEAGE_RDS}'); ids <- names(x); if (is.null(ids) || !length(ids)) ids <- vapply(x, function(i) as.character(i[['lineage_id']]), character(1)); idx <- as.integer('${TASK_ID}'); if (idx < 1 || idx > length(ids)) stop(sprintf('Array task %d is out of range for %d lineages', idx, length(ids)), call. = FALSE); cat(ids[[idx]])")

echo "Resolved lineage_id: ${LINEAGE_ID}"

apptainer exec ${BINDS} "${CONTAINER_URI}" \
  Rscript R/analyze_single_lineage_area.R \
    --lineage_id "${LINEAGE_ID}" \
    --lineage_rds "${LINEAGE_RDS}" \
    --object_sizes_rds "${OBJECT_SIZES_RDS}" \
    --out_dir "${OUT_DIR}"

echo "=========================================================="
echo "Finished lineage ${LINEAGE_ID} at $(date)"
echo "=========================================================="
