import argparse
import os
from pathlib import Path

import numpy as np
from PIL import Image
from tifffile import imread


def normalize_for_display(img: np.ndarray) -> np.ndarray:
    arr = np.asarray(img)
    if arr.ndim == 3:
        arr = arr[..., : min(3, arr.shape[2])].mean(axis=2)
    arr = arr.astype(np.float32, copy=False)
    finite = np.isfinite(arr)
    if not finite.any():
        return np.zeros((arr.shape[0], arr.shape[1]), dtype=np.uint8)

    vals = arr[finite]
    q01, q99 = np.quantile(vals, [0.01, 0.99])
    arr = np.clip(arr, q01, q99)
    arr = arr - np.nanmin(arr)
    vmax = np.nanmax(arr)
    if np.isfinite(vmax) and vmax > 0:
      arr = arr / vmax
    arr = np.clip(arr * 255.0, 0, 255).astype(np.uint8)
    return arr


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert one raw TIFF image to cached PNG.")
    parser.add_argument("--raw-path", required=True)
    parser.add_argument("--out-path", required=True)
    args = parser.parse_args()

    raw_path = Path(args.raw_path)
    out_path = Path(args.out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists():
        print(f"PNG already exists: {out_path}")
        return

    img = imread(raw_path)
    display = normalize_for_display(img)
    Image.fromarray(display).save(out_path)
    print(f"Wrote PNG: {out_path}")


if __name__ == "__main__":
    main()
