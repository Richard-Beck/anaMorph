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
4. Use the dev subset for parser, QC, feature-extraction, and regression tests before running larger jobs on HPC.

## Example usage

Build the manifest:

```bash
Rscript R/index_hpc_images.R \
  --hpc_dir /share/lab_crd/lab_crd/CLONEID/data/LTEEs \
  --passaging_csv core_data/passaging.csv \
  --media_csv core_data/media.csv \
  --out_csv data/image_manifest.csv \
  --out_id_csv data/image_id_summary.csv
```

Select a local dev subset:

```bash
Rscript R/select_dev_subset.R \
  --manifest_csv data/image_manifest.csv \
  --out_image_csv data/dev_subset_images.csv \
  --out_id_csv data/dev_subset_ids.csv \
  --copy_script data/copy_dev_subset.sh \
  --n_groups 4 \
  --ids_per_group 4
```

The subset script writes a bash copy script that can be run on HPC after you adjust the destination path.

## Design intent

- Do not change the canonical HPC folder layout.
- Prefer manifests, summaries, and symlinked working views over file moves.
- Keep only code, metadata, and a small dev subset in the local repo.
- Treat segmentation as an external dependency until the model invocation and outputs are pulled into versioned pipeline code.
