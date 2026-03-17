# anaMorph

This repository is for local development and testing of an image-analysis pipeline for brightfield long-term evolution experiment (LTEE) images. The HPC image store is treated as read-only source data. Local code should index and summarize that data in place rather than reorganizing it.

## Current data contract

- Core experimental metadata lives in [core_data/passaging.csv](./core_data/passaging.csv) and [core_data/media.csv](./core_data/media.csv).
- Raw LTEE images live on the HPC in a single flat directory with filenames like `HGC-27_A10_seedT1_10x_ph_bl.tif`.
- The join key to `passaging.csv$id` is the filename stem before `_10x_ph_<tile>`.
- The tile suffix is one of `bl`, `br`, `tl`, `tr`.
- PNG overlays may also exist with the same stem plus `_overlay`.
- The current external CellposeSAM model is on HPC at `/share/lab_crd/lab_crd/CLONEID/cellpose_segmentation_models/current_model`.

## Recommended workflow

1. Export or refresh database metadata locally with [R/update_db.R](./R/update_db.R) when needed.
2. Run [R/index_hpc_images.R](./R/index_hpc_images.R) against the HPC image directory to build an image manifest joined to passaging/media metadata.
3. Run [R/select_dev_subset.R](./R/select_dev_subset.R) on that manifest to choose a small, reproducible dev subset for local testing and Git sync.
4. Run [scripts/run_dev_subset_cellpose.py](./scripts/run_dev_subset_cellpose.py) on HPC to copy the dev images into the repo and generate tracked segmentation masks.
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

The manifest CSVs, copied dev images, and saved dev masks are intended to be tracked in git and synchronized between local and HPC copies of the repo. Regenerable summary tables are not.

Run CellposeSAM on the dev subset:

```bash
python scripts/run_dev_subset_cellpose.py \
  --manifest-csv manifests/dev_subset_images.csv \
  --pretrained-model /share/lab_crd/lab_crd/CLONEID/cellpose_segmentation_models/current_model \
  --raw-dir dev_data/raw \
  --mask-dir dev_data/segmentation_masks \
  --run-manifest manifests/dev_subset_segmentation_runs.csv \
  --image-manifest manifests/dev_subset_segmentation_images.csv
```

This uses only:

- the dev subset images
- `gpu=True`
- `pretrained_model=/share/lab_crd/lab_crd/CLONEID/cellpose_segmentation_models/current_model`
- `model.eval(imgs)`

Generate summaries later when needed:

```bash
python scripts/summarize_dev_subset_masks.py \
  --segmentation-manifest manifests/dev_subset_segmentation_images.csv \
  --image-summary-csv data/dev_subset_image_summary.csv \
  --object-summary-csv data/dev_subset_object_summary.csv
```

## Design intent

- Do not change the canonical HPC folder layout.
- Prefer manifests, summaries, and symlinked working views over file moves.
- Keep only code, metadata, the small dev subset, and small segmentation artifacts in the local repo.
- Keep syncable manifests and subset definitions under tracked paths in the repo.
- Treat segmentation as an external dependency until the model invocation and outputs are pulled into versioned pipeline code.
- Keep large regenerable object summaries out of version control when possible.
