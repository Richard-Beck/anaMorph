import argparse
import os

import numpy as np
import pandas as pd
from skimage.measure import regionprops_table
from tifffile import imread


def parse_args():
    parser = argparse.ArgumentParser(description="Generate image and object summaries from tracked dev images and masks.")
    parser.add_argument("--segmentation-manifest", default="manifests/dev_subset_segmentation_images.csv")
    parser.add_argument("--image-summary-csv", default="data/dev_subset_image_summary.csv")
    parser.add_argument("--object-summary-csv", default="data/dev_subset_object_summary.csv")
    return parser.parse_args()


def ensure_parent(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def normalize_relpath(path):
    return path.replace("\\", "/")


def main():
    args = parse_args()
    manifest = pd.read_csv(args.segmentation_manifest)
    if manifest.empty:
        raise SystemExit("No rows found in segmentation manifest.")

    image_rows = []
    object_rows = []

    for row in manifest.to_dict(orient="records"):
        image = imread(row["raw_relpath"])
        mask = np.load(row["mask_relpath"])["mask"]

        object_count = int(np.max(mask)) if np.size(mask) else 0
        foreground_pixels = int(np.sum(mask > 0))
        image_rows.append(
            {
                "run_id": row["run_id"],
                "id": row.get("id", ""),
                "field": row.get("field", ""),
                "filename": row["filename"],
                "raw_relpath": normalize_relpath(row["raw_relpath"]),
                "mask_relpath": normalize_relpath(row["mask_relpath"]),
                "height_px": int(mask.shape[0]),
                "width_px": int(mask.shape[1]),
                "object_count": object_count,
                "foreground_fraction": foreground_pixels / mask.size if mask.size else 0.0,
                "intensity_mean": float(np.mean(image)),
                "intensity_sd": float(np.std(image)),
            }
        )

        if object_count == 0:
            continue

        props = regionprops_table(
            mask,
            intensity_image=image,
            properties=(
                "label",
                "area",
                "bbox",
                "centroid",
                "axis_major_length",
                "axis_minor_length",
                "eccentricity",
                "mean_intensity",
                "intensity_max",
                "intensity_min",
            ),
        )
        props_df = pd.DataFrame(props)
        if props_df.empty:
            continue

        props_df["aspect_ratio"] = props_df["axis_major_length"] / props_df["axis_minor_length"].replace(0, np.nan)
        props_df["run_id"] = row["run_id"]
        props_df["id"] = row.get("id", "")
        props_df["field"] = row.get("field", "")
        props_df["filename"] = row["filename"]
        props_df["raw_relpath"] = normalize_relpath(row["raw_relpath"])
        props_df["mask_relpath"] = normalize_relpath(row["mask_relpath"])
        object_rows.append(props_df)

    ensure_parent(args.image_summary_csv)
    ensure_parent(args.object_summary_csv)

    pd.DataFrame(image_rows).to_csv(args.image_summary_csv, index=False)
    if object_rows:
        pd.concat(object_rows, ignore_index=True).to_csv(args.object_summary_csv, index=False)
    else:
        pd.DataFrame().to_csv(args.object_summary_csv, index=False)

    print(f"Wrote image summary: {args.image_summary_csv}")
    print(f"Wrote object summary: {args.object_summary_csv}")


if __name__ == "__main__":
    main()
