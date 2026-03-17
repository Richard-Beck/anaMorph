import argparse
import csv
import os
import shutil
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from cellpose import io, models
from cellpose.io import imread


def parse_args():
    parser = argparse.ArgumentParser(description="Run CellposeSAM on the tracked dev subset.")
    parser.add_argument("--manifest-csv", default="manifests/dev_subset_images.csv")
    parser.add_argument(
        "--pretrained-model",
        default="/share/lab_crd/lab_crd/CLONEID/cellpose_segmentation_models/current_model",
    )
    parser.add_argument("--raw-dir", default="dev_data/raw")
    parser.add_argument("--mask-dir", default="dev_data/segmentation_masks")
    parser.add_argument("--run-manifest", default="manifests/dev_subset_segmentation_runs.csv")
    parser.add_argument("--image-manifest", default="manifests/dev_subset_segmentation_images.csv")
    parser.add_argument("--run-label", default="")
    return parser.parse_args()


def ensure_parent(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def normalize_relpath(path):
    return path.replace("\\", "/")


def main():
    args = parse_args()
    io.logger_setup()

    manifest = pd.read_csv(args.manifest_csv)
    if manifest.empty:
        raise SystemExit("No rows found in manifest.")

    files = manifest["filepath"].tolist()
    imgs = [imread(path) for path in files]

    model = models.CellposeModel(gpu=True, pretrained_model=args.pretrained_model)
    masks, flows, styles = model.eval(imgs)

    os.makedirs(args.raw_dir, exist_ok=True)
    os.makedirs(args.mask_dir, exist_ok=True)
    ensure_parent(args.run_manifest)
    ensure_parent(args.image_manifest)

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    run_id = args.run_label or f"dev_subset_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"

    image_rows = []
    for row, image, mask in zip(manifest.to_dict(orient="records"), imgs, masks):
        source_path = row["filepath"]
        source_name = os.path.basename(source_path)
        stem, _ = os.path.splitext(source_name)

        raw_dest = os.path.join(args.raw_dir, source_name)
        mask_dest = os.path.join(args.mask_dir, f"{stem}_mask.npz")

        if not os.path.exists(raw_dest):
            shutil.copy2(source_path, raw_dest)

        np.savez_compressed(mask_dest, mask=np.asarray(mask, dtype=np.int32))

        image_rows.append(
            {
                "run_id": run_id,
                "timestamp_utc": timestamp,
                "id": row.get("id", ""),
                "field": row.get("field", ""),
                "filename": source_name,
                "source_filepath": normalize_relpath(source_path),
                "raw_relpath": normalize_relpath(raw_dest),
                "mask_relpath": normalize_relpath(mask_dest),
                "image_shape": "x".join(str(x) for x in np.asarray(image).shape),
                "mask_shape": "x".join(str(x) for x in np.asarray(mask).shape),
                "object_count": int(np.max(mask)) if np.size(mask) else 0,
            }
        )

    run_row = {
        "run_id": run_id,
        "timestamp_utc": timestamp,
        "manifest_csv": normalize_relpath(args.manifest_csv),
        "pretrained_model": args.pretrained_model,
        "gpu": True,
        "image_count": len(image_rows),
        "raw_dir": normalize_relpath(args.raw_dir),
        "mask_dir": normalize_relpath(args.mask_dir),
    }

    with open(args.run_manifest, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(run_row.keys()))
        writer.writeheader()
        writer.writerow(run_row)

    image_df = pd.DataFrame(image_rows)
    image_df.to_csv(args.image_manifest, index=False)

    print(f"Run ID: {run_id}")
    print(f"Segmented images: {len(image_rows)}")
    print(f"Wrote run manifest: {args.run_manifest}")
    print(f"Wrote image manifest: {args.image_manifest}")


if __name__ == "__main__":
    main()
