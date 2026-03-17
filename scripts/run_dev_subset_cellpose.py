import argparse
import csv
import os
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from cellpose import io, models
from cellpose.io import imread
from PIL import Image
from skimage.measure import regionprops
from tifffile import imwrite
from torchvision import models as tv_models


def parse_args():
    parser = argparse.ArgumentParser(description="Run CellposeSAM on the tracked dev subset.")
    parser.add_argument("--manifest-csv", default="manifests/dev_subset_images.csv")
    parser.add_argument(
        "--pretrained-model",
        default="/share/lab_crd/lab_crd/CLONEID/cellpose_segmentation_models/current_model",
    )
    parser.add_argument("--raw-dir", default="dev_data/raw")
    parser.add_argument("--mask-dir", default="dev_data/segmentation_masks")
    parser.add_argument("--embedding-dir", default="dev_data/embeddings")
    parser.add_argument("--run-manifest", default="manifests/dev_subset_segmentation_runs.csv")
    parser.add_argument("--image-manifest", default="manifests/dev_subset_segmentation_images.csv")
    parser.add_argument("--run-label", default="")
    parser.add_argument("--preferred-field", default="bl")
    parser.add_argument("--crop-height", type=int, default=1100)
    parser.add_argument("--crop-width", type=int, default=1400)
    return parser.parse_args()


def ensure_parent(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def normalize_relpath(path):
    return path.replace("\\", "/")


def mask_output_path(mask_dir, source_name):
    stem, _ = os.path.splitext(source_name)
    return os.path.join(mask_dir, f"{stem}_mask.tif")


def embedding_output_path(embedding_dir, source_name):
    stem, _ = os.path.splitext(source_name)
    return os.path.join(embedding_dir, f"{stem}_embeddings.csv")


def crop_top_left(image, crop_height, crop_width):
    y_max = min(crop_height, image.shape[0])
    x_max = min(crop_width, image.shape[1])
    if image.ndim == 2:
        return image[:y_max, :x_max]
    return image[:y_max, :x_max, ...]


def select_one_field_per_id(manifest, preferred_field):
    field_order = [preferred_field, "bl", "br", "tl", "tr"]
    priority = {field: idx for idx, field in enumerate(field_order)}

    manifest = manifest.copy()
    manifest["field_priority"] = manifest["field"].map(lambda x: priority.get(x, len(priority)))
    manifest = manifest.sort_values(["id", "field_priority", "filename"]).drop_duplicates(subset=["id"], keep="first")
    return manifest.drop(columns=["field_priority"])


def build_embedding_model():
    weights = tv_models.ResNet18_Weights.DEFAULT
    model = tv_models.resnet18(weights=weights)
    model.fc = nn.Identity()
    model.eval()
    return model, weights.transforms()


def to_grayscale_uint8(image):
    if image.ndim == 2:
        if image.dtype == np.uint8:
            return image
        return np.asarray(Image.fromarray(image).convert("L"))

    pil_image = Image.fromarray(image)
    return np.asarray(pil_image.convert("L"))


def object_embedding_rows(image, mask, model, preprocess, device):
    gray_image = to_grayscale_uint8(image)
    embedding_columns = [f"embedding_{idx:03d}" for idx in range(512)]
    base_columns = [
        "label",
        "area_px",
        "log1p_area",
        "bbox_min_row",
        "bbox_min_col",
        "bbox_max_row",
        "bbox_max_col",
    ]
    rows = []

    for prop in regionprops(mask):
        min_row, min_col, max_row, max_col = prop.bbox
        gray_crop = gray_image[min_row:max_row, min_col:max_col].copy()
        gray_crop[~prop.image] = 0

        pil_crop = Image.fromarray(gray_crop, mode="L").convert("RGB")
        x = preprocess(pil_crop).unsqueeze(0).to(device)

        with torch.inference_mode():
            emb = model(x).squeeze(0).cpu().numpy()

        row = {
            "label": int(prop.label),
            "area_px": int(prop.area),
            "log1p_area": float(np.log1p(prop.area)),
            "bbox_min_row": int(min_row),
            "bbox_min_col": int(min_col),
            "bbox_max_row": int(max_row),
            "bbox_max_col": int(max_col),
        }
        row.update({col: float(val) for col, val in zip(embedding_columns, emb)})
        rows.append(row)

    return rows, base_columns + embedding_columns


def main():
    args = parse_args()
    io.logger_setup()

    manifest = pd.read_csv(args.manifest_csv)
    if manifest.empty:
        raise SystemExit("No rows found in manifest.")
    manifest = select_one_field_per_id(manifest, args.preferred_field)

    files = manifest["filepath"].tolist()
    source_imgs = [imread(path) for path in files]
    imgs = [crop_top_left(img, args.crop_height, args.crop_width) for img in source_imgs]

    model = models.CellposeModel(gpu=True, pretrained_model=args.pretrained_model)
    masks, flows, styles = model.eval(imgs)
    embedding_model, preprocess = build_embedding_model()
    embedding_device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    embedding_model = embedding_model.to(embedding_device)

    os.makedirs(args.raw_dir, exist_ok=True)
    os.makedirs(args.mask_dir, exist_ok=True)
    os.makedirs(args.embedding_dir, exist_ok=True)
    ensure_parent(args.run_manifest)
    ensure_parent(args.image_manifest)

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    run_id = args.run_label or f"dev_subset_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"

    image_rows = []
    for row, source_image, image, mask in zip(manifest.to_dict(orient="records"), source_imgs, imgs, masks):
        source_path = row["filepath"]
        source_name = os.path.basename(source_path)

        raw_dest = os.path.join(args.raw_dir, source_name)
        mask_dest = mask_output_path(args.mask_dir, source_name)
        embedding_dest = embedding_output_path(args.embedding_dir, source_name)

        io.imsave(raw_dest, image)

        imwrite(mask_dest, np.asarray(mask, dtype=np.int32))
        embedding_rows, embedding_columns = object_embedding_rows(
            image=image,
            mask=np.asarray(mask, dtype=np.int32),
            model=embedding_model,
            preprocess=preprocess,
            device=embedding_device,
        )
        pd.DataFrame(embedding_rows, columns=embedding_columns).to_csv(embedding_dest, index=False)

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
                "embedding_relpath": normalize_relpath(embedding_dest),
                "source_image_shape": "x".join(str(x) for x in np.asarray(source_image).shape),
                "image_shape": "x".join(str(x) for x in np.asarray(image).shape),
                "mask_shape": "x".join(str(x) for x in np.asarray(mask).shape),
                "object_count": int(np.max(mask)) if np.size(mask) else 0,
                "embedding_count": len(embedding_rows),
                "embedding_dim": 512,
                "crop_height": args.crop_height,
                "crop_width": args.crop_width,
            }
        )

    run_row = {
        "run_id": run_id,
        "timestamp_utc": timestamp,
        "manifest_csv": normalize_relpath(args.manifest_csv),
        "pretrained_model": args.pretrained_model,
        "gpu": True,
        "preferred_field": args.preferred_field,
        "image_count": len(image_rows),
        "raw_dir": normalize_relpath(args.raw_dir),
        "mask_dir": normalize_relpath(args.mask_dir),
        "embedding_dir": normalize_relpath(args.embedding_dir),
        "embedding_model": "resnet18",
        "embedding_weights": "ResNet18_Weights.DEFAULT",
        "embedding_dim": 512,
        "crop_height": args.crop_height,
        "crop_width": args.crop_width,
        "crop_origin": "top_left",
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
