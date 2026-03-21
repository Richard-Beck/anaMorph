#!/bin/bash

#----------------------------------------------------------------#
# SBATCH - SLURM Resource Directives
#----------------------------------------------------------------#
#SBATCH --job-name=raw_img_embed
#SBATCH --output=logs/raw-img-embed-%j.out
#SBATCH --error=logs/raw-img-embed-%j.err

#SBATCH --time=24:00:00
#SBATCH --qos=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --gres=gpu:a30:1

#----------------------------------------------------------------#
# Environment Setup & Job Execution
#----------------------------------------------------------------#

set -eo pipefail

echo "=========================================================="
echo "Job Started at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURMD_NODENAME"
echo "=========================================================="

mkdir -p logs
mkdir -p data/all_images/manifests

source /app/eb/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh
set +u
conda activate cpose
set -u

echo "Python environment activated. Running the script..."
python scripts/run_all_raw_image_embeddings.py \
  --raw-data-dir /share/lab_crd/lab_crd/CLONEID/data/LTEEs \
  --output-csv data/all_images/manifests/raw_image_embeddings.csv \
  --run-manifest data/all_images/manifests/raw_image_embedding_runs.csv \
  --embedding-batch-size 64

echo "=========================================================="
echo "Job Finished at $(date)"
echo "=========================================================="
