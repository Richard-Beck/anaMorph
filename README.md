# anaMorph

This repository is for local development and testing of an image-analysis pipeline for brightfield long-term evolution experiment (LTEE) images. The HPC image store is treated as read-only source data. Local code should index and summarize that data in place rather than reorganizing it.

## Current data contract

- Core experimental metadata lives in [core_data/passaging.csv](./core_data/passaging.csv) and [core_data/media.csv](./core_data/media.csv).
- Raw LTEE images live on the HPC in a single flat directory with filenames like `HGC-27_A10_seedT1_10x_ph_bl.tif`.
- The join key to `passaging.csv$id` is the filename stem before `_10x_ph_<tile>`.
- The tile suffix is one of `bl`, `br`, `tl`, `tr`.
- PNG overlays may also exist with the same stem plus `_overlay`.
- The current external CellposeSAM model is on HPC at `/share/lab_crd/lab_crd/CLONEID/cellpose_segmentation_models/current_model`.
- Current working assumption: a fixed top-left crop retaining `y < 1100` and `x < 1400` removes the scale bar across the dev subset, so tracked dev raws and downstream segmentations should be generated from that crop rather than the full frame.

## Recommended workflow

1. Export or refresh database metadata locally with [R/update_db.R](./R/update_db.R) when needed.
2. Run [R/index_hpc_images.R](./R/index_hpc_images.R) against the HPC image directory to build an image manifest joined to passaging/media metadata.
3. Run [R/select_dev_subset.R](./R/select_dev_subset.R) on that manifest to choose a small, reproducible dev subset for local testing and Git sync.
4. Run [scripts/run_dev_subset_cellpose.py](./scripts/run_dev_subset_cellpose.py) on HPC to copy the dev images into the repo and generate tracked segmentation masks as labeled TIFFs.
   By default it segments one image per `id` and prefers the `bl` field when multiple fields exist.
   It also crops each source image to the top-left `1100 x 1400` region before saving the tracked dev raw image and running Cellpose, to exclude the scale bar artifact.
   For each segmented image, it exports one CSV in `dev_data/embeddings/` with one row per object containing the object label, bbox, area, `log1p(area)`, and a 512-d embedding from pretrained `ResNet18_Weights.DEFAULT`.
   The current wrapper streams work in small batches rather than loading the whole dev subset at once: Cellpose runs on image batches and embeddings run on per-image object batches, then outputs are written immediately.
5. Generate image/object summaries on demand with [scripts/summarize_dev_subset_masks.py](./scripts/summarize_dev_subset_masks.py) instead of committing large object tables.

## Example usage

Build the manifest:

```bash
Rscript R/index_hpc_images.R \
  --hpc_dir /share/lab_crd/lab_crd/CLONEID/data/LTEEs \
  --passaging_csv core_data/passaging.csv \
  --media_csv core_data/media.csv \
  --out_csv manifests/image_manifest.csv \
  --out_id_csv manifests/image_id_summary.csv
```

Select a local dev subset:

```bash
Rscript R/select_dev_subset.R \
  --manifest_csv manifests/image_manifest.csv \
  --out_image_csv manifests/dev_subset_images.csv \
  --out_id_csv manifests/dev_subset_ids.csv \
  --copy_script scripts/copy_dev_subset.sh \
  --n_groups 4 \
  --ids_per_group 4
```

The manifest CSVs, copied dev images, and saved TIFF masks are intended to be tracked in git and synchronized between local and HPC copies of the repo. Regenerable summary tables are not.

Run CellposeSAM on the dev subset:

```bash
python scripts/run_dev_subset_cellpose.py \
  --manifest-csv manifests/dev_subset_images.csv \
  --pretrained-model /share/lab_crd/lab_crd/CLONEID/cellpose_segmentation_models/current_model \
  --raw-dir dev_data/raw \
  --mask-dir dev_data/segmentation_masks \
  --embedding-dir dev_data/embeddings \
  --run-manifest manifests/dev_subset_segmentation_runs.csv \
  --image-manifest manifests/dev_subset_segmentation_images.csv \
  --preferred-field bl \
  --cellpose-batch-size 8 \
  --embedding-batch-size 256
```

This uses only:

- the dev subset images
- `gpu=True`
- `pretrained_model=/share/lab_crd/lab_crd/CLONEID/cellpose_segmentation_models/current_model`
- `model.eval(imgs)`
- pretrained `ResNet18_Weights.DEFAULT` embeddings for each segmented object, saved per image as CSV

Generate summaries later when needed:

```bash
python scripts/summarize_dev_subset_masks.py \
  --segmentation-manifest manifests/dev_subset_segmentation_images.csv \
  --image-summary-csv data/dev_subset_image_summary.csv \
  --object-summary-csv data/dev_subset_object_summary.csv
```

Reusable embedding clustering helpers live in [R/embedding_clustering_utils.R](./R/embedding_clustering_utils.R). Example:

```r
source("R/embedding_clustering_utils.R")

image_names <- readr::read_csv("manifests/dev_subset_images.csv", show_col_types = FALSE)$filename
embedding_tbl <- load_embedding_vectors(image_names)
pca_result <- run_embedding_pca(embedding_tbl, n_pcs = 20)
cluster_result <- cluster_reduced_embeddings(pca_result$scores, feature_cols = pca_result$score_cols)
```

## Full HPC Run

For full-dataset processing on HPC, use [scripts/run_all_raw_cellpose.py](./scripts/run_all_raw_cellpose.py). It scans the raw LTEE image directory directly, applies the same top-left crop (`y < 1100`, `x < 1400`), writes segmentation masks and per-image embedding CSVs under `data/all_images/`, and does not save cropped raw images.

Example direct run:

```bash
python scripts/run_all_raw_cellpose.py \
  --raw-data-dir /share/lab_crd/lab_crd/CLONEID/data/LTEEs \
  --pretrained-model /share/lab_crd/lab_crd/CLONEID/cellpose_segmentation_models/current_model \
  --mask-dir data/all_images/segmentation_masks \
  --embedding-dir data/all_images/embeddings \
  --run-manifest data/all_images/manifests/segmentation_runs.csv \
  --image-manifest data/all_images/manifests/segmentation_images.csv \
  --cellpose-batch-size 16 \
  --embedding-batch-size 512
```

For cluster submission, use [scripts/submit_all_raw_cellpose.slurm](./scripts/submit_all_raw_cellpose.slurm).

To save a compact per-image list of object sizes from the full-run embedding CSVs:

```bash
Rscript R/build_all_image_object_sizes.R \
  --embedding_dir data/all_images/embeddings \
  --out_rds data/all_images/object_sizes.rds
```

To analyze cell-area residual distributions for one lineage at a time:

```bash
Rscript R/analyze_single_lineage_area.R \
  --lineage_id 7 \
  --lineage_rds core_data/lineages.Rds \
  --object_sizes_rds data/all_images/object_sizes.rds \
  --out_dir data/lineage_area
```

This writes one lineage-specific folder with per-cell residuals plus image- and passage-level summaries after a small-object mixture adjustment and growth-state normalization.

To submit the single-lineage analysis across all saved lineages on SLURM as an array job:

```bash
sbatch scripts/submit_lineage_area_analysis.slurm \
  core_data/lineages.Rds \
  data/all_images/object_sizes.rds \
  data/lineage_area
```

The wrapper resolves one `lineage_id` per array task from `core_data/lineages.Rds` and runs `R/analyze_single_lineage_area.R`. The current default array range is `1-88`; if the lineage count changes, update the `#SBATCH --array=` line or override it at submission time with `sbatch --array=1-N ...`.

## Design intent

- Do not change the canonical HPC folder layout.
- Prefer manifests, summaries, and symlinked working views over file moves.
- Keep only code, metadata, the small dev subset, and small segmentation artifacts in the local repo.
- Keep syncable manifests and subset definitions under tracked paths in the repo.
- Treat segmentation as an external dependency until the model invocation and outputs are pulled into versioned pipeline code.
- Keep large regenerable object summaries out of version control when possible.

## Dev Notes

- Current performance issue: adding per-object ResNet-18 embeddings makes `scripts/run_dev_subset_cellpose.py` unacceptably slow for routine iteration.
- Observed behavior: GPU utilization remains low during the embedding stage, so the current implementation is not using hardware efficiently.
- Pipeline constraint: doing Cellpose for all images first and deferring all embedding work until afterward is not viable at larger scale because the intermediate state will exceed available memory.
- Working assumption for future optimization: the main bottleneck is likely per-object serial preprocessing/inference rather than Cellpose segmentation itself, so batching and a more GPU-efficient embedding path should be prioritized.
