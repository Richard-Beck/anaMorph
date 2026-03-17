#!/bin/bash

#----------------------------------------------------------------#
# SBATCH - SLURM Resource Directives
#----------------------------------------------------------------#
#SBATCH --job-name=all_raw_cellpose
#SBATCH --output=logs/all-raw-cellpose-%j.out
#SBATCH --error=logs/all-raw-cellpose-%j.err

#SBATCH --time=24:00:00
#SBATCH --qos=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=128G
#SBATCH --gres=gpu:a30:1

#----------------------------------------------------------------#
# Environment Setup & Job Execution
#----------------------------------------------------------------#

set -euo pipefail

echo "=========================================================="
echo "Job Started at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURMD_NODENAME"
echo "=========================================================="

mkdir -p logs
mkdir -p data/all_images/manifests

source /app/eb/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh
conda activate cpose

echo "Python environment activated. Running the script..."
python scripts/run_all_raw_cellpose.py \
  --raw-data-dir /share/lab_crd/lab_crd/CLONEID/data/LTEEs \
  --pretrained-model /share/lab_crd/lab_crd/CLONEID/cellpose_segmentation_models/current_model \
  --mask-dir data/all_images/segmentation_masks \
  --embedding-dir data/all_images/embeddings \
  --run-manifest data/all_images/manifests/segmentation_runs.csv \
  --image-manifest data/all_images/manifests/segmentation_images.csv \
  --cellpose-batch-size 16 \
  --embedding-batch-size 512

echo "=========================================================="
echo "Job Finished at $(date)"
echo "=========================================================="
