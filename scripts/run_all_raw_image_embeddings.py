import argparse
import csv
import os
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from cellpose.io import imread
from PIL import Image
from torchvision import models as tv_models


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate one ResNet18 embedding vector for each raw LTEE image."
    )
    parser.add_argument(
        "--raw-data-dir",
        default="/share/lab_crd/lab_crd/CLONEID/data/LTEEs",
        help="Directory containing raw TIFF images.",
    )
    parser.add_argument(
        "--output-csv",
        default="data/all_images/manifests/raw_image_embeddings.csv",
        help="Single CSV that stores one embedding row per raw image.",
    )
    parser.add_argument(
        "--run-manifest",
        default="data/all_images/manifests/raw_image_embedding_runs.csv",
        help="CSV log with one row per embedding run.",
    )
    parser.add_argument("--run-label", default="")
    parser.add_argument("--crop-height", type=int, default=1100)
    parser.add_argument("--crop-width", type=int, default=1400)
    parser.add_argument("--embedding-batch-size", type=int, default=64)
    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        help="Rebuild the output CSV from scratch instead of appending only missing images.",
    )
    return parser.parse_args()


def ensure_parent(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def normalize_relpath(path):
    return path.replace("\\", "/")


def iter_batches(items, batch_size):
    for start in range(0, len(items), batch_size):
        yield items[start : start + batch_size]


def list_raw_images(raw_data_dir):
    exts = {".tif", ".tiff"}
    filepaths = []
    for name in os.listdir(raw_data_dir):
        full_path = os.path.join(raw_data_dir, name)
        if not os.path.isfile(full_path):
            continue
        stem, ext = os.path.splitext(name)
        if ext.lower() not in exts:
            continue
        if stem.endswith("_overlay"):
            continue
        filepaths.append(full_path)
    filepaths.sort()
    return filepaths


def parse_image_metadata(filepath):
    filename = os.path.basename(filepath)
    stem, _ = os.path.splitext(filename)
    field = ""
    image_id = stem
    marker = "_10x_ph_"
    if marker in stem:
        image_id, field = stem.rsplit(marker, 1)
    return {
        "id": image_id,
        "field": field,
        "filename": filename,
        "filepath": filepath,
    }


def crop_top_left(image, crop_height, crop_width):
    y_max = min(crop_height, image.shape[0])
    x_max = min(crop_width, image.shape[1])
    if image.ndim == 2:
        return image[:y_max, :x_max]
    return image[:y_max, :x_max, ...]


def build_embedding_model():
    weights = tv_models.ResNet18_Weights.DEFAULT
    model = tv_models.resnet18(weights=weights)
    model.fc = nn.Identity()
    model.eval()
    return model, weights.transforms()


def coerce_to_uint8(image):
    arr = np.asarray(image)
    if arr.dtype == np.uint8:
        return arr
    if arr.dtype == np.bool_:
        return arr.astype(np.uint8) * 255

    arr = arr.astype(np.float32, copy=False)
    finite_mask = np.isfinite(arr)
    if not finite_mask.any():
        return np.zeros(arr.shape, dtype=np.uint8)

    finite_vals = arr[finite_mask]
    min_val = float(finite_vals.min())
    max_val = float(finite_vals.max())
    if max_val <= min_val:
        out = np.zeros(arr.shape, dtype=np.uint8)
        out[finite_mask] = int(np.clip(min_val, 0, 255))
        return out

    scaled = (arr - min_val) / (max_val - min_val)
    scaled = np.clip(scaled * 255.0, 0, 255)
    scaled[~finite_mask] = 0
    return scaled.astype(np.uint8)


def prepare_rgb_pil(image):
    arr = coerce_to_uint8(image)
    if arr.ndim == 2:
        return Image.fromarray(arr, mode="L").convert("RGB")
    if arr.ndim == 3 and arr.shape[2] == 1:
        return Image.fromarray(arr[:, :, 0], mode="L").convert("RGB")
    if arr.ndim == 3 and arr.shape[2] == 2:
        padded = np.concatenate([arr, arr[:, :, :1]], axis=2)
        return Image.fromarray(padded, mode="RGB")
    if arr.ndim == 3 and arr.shape[2] >= 3:
        return Image.fromarray(arr[:, :, :3], mode="RGB")
    raise ValueError(f"Unsupported image shape for embedding: {arr.shape}")


def append_rows(csv_path, rows):
    if not rows:
        return
    ensure_parent(csv_path)
    file_exists = os.path.exists(csv_path) and os.path.getsize(csv_path) > 0
    with open(csv_path, "a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        if not file_exists:
            writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def load_processed_sources(output_csv):
    if not os.path.exists(output_csv):
        return set()
    try:
        existing = pd.read_csv(output_csv, usecols=["source_filepath"])
    except (ValueError, pd.errors.EmptyDataError):
        return set()
    if existing.empty:
        return set()
    return set(existing["source_filepath"].dropna().astype(str))


def remove_existing_output(output_csv):
    if os.path.exists(output_csv):
        os.remove(output_csv)


def embedding_rows_for_batch(
    batch_rows,
    cropped_images,
    source_images,
    model,
    preprocess,
    device,
    run_id,
    timestamp,
    crop_height,
    crop_width,
):
    embedding_columns = [f"embedding_{idx:03d}" for idx in range(512)]
    tensors = [preprocess(prepare_rgb_pil(image)) for image in cropped_images]
    x = torch.stack(tensors, dim=0).to(device, non_blocking=True)
    with torch.inference_mode():
        batch_embeddings = model(x).cpu().numpy()

    output_rows = []
    for row, source_image, cropped_image, embedding in zip(batch_rows, source_images, cropped_images, batch_embeddings):
        output_row = {
            "id": row["id"],
            "field": row["field"],
            "filename": row["filename"],
            "source_filepath": normalize_relpath(row["filepath"]),
            "source_image_shape": "x".join(str(x) for x in np.asarray(source_image).shape),
            "image_shape": "x".join(str(x) for x in np.asarray(cropped_image).shape),
            "embedding_dim": 512,
            "crop_height": crop_height,
            "crop_width": crop_width,
            "run_id": run_id,
            "timestamp_utc": timestamp,
        }
        output_row.update({col: float(val) for col, val in zip(embedding_columns, embedding)})
        output_rows.append(output_row)
    return output_rows


def append_manifest_row(csv_path, row):
    append_rows(csv_path, [row])


def main():
    args = parse_args()

    filepaths = list_raw_images(args.raw_data_dir)
    if not filepaths:
        raise SystemExit(f"No raw TIFF images found in {args.raw_data_dir}")

    if args.overwrite_existing:
        remove_existing_output(args.output_csv)
        processed_sources = set()
    else:
        processed_sources = load_processed_sources(args.output_csv)

    manifest_rows = [parse_image_metadata(path) for path in filepaths]
    manifest_rows = [
        row for row in manifest_rows if normalize_relpath(row["filepath"]) not in processed_sources
    ]

    if not manifest_rows:
        print("No new raw images to embed.")
        print(f"Existing output CSV: {args.output_csv}")
        return

    embedding_model, preprocess = build_embedding_model()
    embedding_device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    embedding_model = embedding_model.to(embedding_device)
    embedding_model.eval()

    ensure_parent(args.output_csv)
    ensure_parent(args.run_manifest)

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    run_id = args.run_label or f"raw_image_embeddings_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
    processed_count = 0

    for manifest_batch in iter_batches(manifest_rows, args.embedding_batch_size):
        source_images = [imread(row["filepath"]) for row in manifest_batch]
        cropped_images = [crop_top_left(img, args.crop_height, args.crop_width) for img in source_images]
        output_rows = embedding_rows_for_batch(
            batch_rows=manifest_batch,
            cropped_images=cropped_images,
            source_images=source_images,
            model=embedding_model,
            preprocess=preprocess,
            device=embedding_device,
            run_id=run_id,
            timestamp=timestamp,
            crop_height=args.crop_height,
            crop_width=args.crop_width,
        )
        append_rows(args.output_csv, output_rows)
        processed_count += len(output_rows)
        print(f"Appended {len(output_rows)} image embeddings. Total written this run: {processed_count}")

    run_row = {
        "run_id": run_id,
        "timestamp_utc": timestamp,
        "raw_data_dir": normalize_relpath(args.raw_data_dir),
        "output_csv": normalize_relpath(args.output_csv),
        "image_count": processed_count,
        "embedding_model": "resnet18",
        "embedding_weights": "ResNet18_Weights.DEFAULT",
        "embedding_dim": 512,
        "embedding_batch_size": args.embedding_batch_size,
        "crop_height": args.crop_height,
        "crop_width": args.crop_width,
        "crop_origin": "top_left",
        "device": str(embedding_device),
        "overwrite_existing": bool(args.overwrite_existing),
    }
    append_manifest_row(args.run_manifest, run_row)

    print(f"Run ID: {run_id}")
    print(f"Embedded images: {processed_count}")
    print(f"Wrote output CSV: {args.output_csv}")
    print(f"Appended run manifest: {args.run_manifest}")


if __name__ == "__main__":
    main()
