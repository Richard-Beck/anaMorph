#!/usr/bin/env python3
import argparse
import base64
import csv
import html
from pathlib import Path
from typing import Optional


PALETTE = [
    "#d60000",
    "#018700",
    "#b500ff",
    "#05acc6",
    "#97ff00",
    "#ffa52f",
    "#ff8ec8",
    "#79525e",
    "#00fdcf",
    "#afa5ff",
    "#93ac83",
    "#9a6900",
    "#366962",
    "#d18fff",
    "#fdf490",
    "#c86e66",
    "#9ee2ff",
    "#00c846",
    "#f4b2a9",
    "#9c0000",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a lineage stack HTML with cluster-colored bounding box overlays."
    )
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--lineage-id", required=True)
    parser.add_argument("--field", default="tl")
    parser.add_argument(
        "--image-summary",
        default=None,
        help="Optional explicit image_summary.csv path. Defaults to data/lineage_area/lineage_<id>/image_summary.csv",
    )
    parser.add_argument(
        "--cluster-csv",
        default=None,
        help="Optional explicit best_clustering.csv path. Defaults to pooled clustering output for the lineage.",
    )
    parser.add_argument(
        "--segmentation-manifest",
        default="data/all_images/manifests/segmentation_images.csv",
        help="Segmentation manifest with source and cropped image shapes.",
    )
    parser.add_argument(
        "--embedding-dir",
        default="data/all_images/embeddings",
        help="Embedding CSV directory, relative to project root unless absolute.",
    )
    parser.add_argument(
        "--png-dir",
        default=None,
        help="Optional explicit lineage PNG directory. Defaults to data/lineage_png/lineage_<id>/",
    )
    parser.add_argument(
        "--out-html",
        required=True,
        help="Output HTML path, relative to project root unless absolute.",
    )
    return parser.parse_args()


def resolve_path(project_root: Path, value: Optional[str]) -> Optional[Path]:
    if value is None:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return project_root / path


def lineage_dir_name(lineage_id: str) -> str:
    return f"lineage_{lineage_id}"


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_int(value: str) -> int:
    return int(float(value))


def as_float(value: str) -> float:
    return float(value)


def load_stack_rows(
    image_summary_path: Path,
    png_dir: Path,
    field: str,
    image_dims: dict[str, dict[str, int]],
) -> list[dict[str, object]]:
    rows = []
    for row in read_csv_rows(image_summary_path):
        if row.get("field") != field:
            continue
        filename = row["filename"]
        dims = image_dims.get(filename)
        if dims is None:
            continue
        png_name = Path(filename).with_suffix(".png").name
        png_path = png_dir / png_name
        if not png_path.exists():
            continue
        rows.append(
            {
                "filename": filename,
                "png_name": png_name,
                "png_path": png_path,
                "passage_position": as_int(row["passage_position"]),
                "passage_id": row["passage_id"],
                "time_hours": as_float(row["mean_time_since_passage_start_hours"]),
                **dims,
            }
        )
    rows.sort(key=lambda row: (row["passage_position"], row["time_hours"], row["png_name"]))
    return rows


def load_cluster_boxes(
    cluster_csv: Path,
    embedding_dir: Path,
    keep_filenames: set[str],
) -> dict[str, list[dict[str, object]]]:
    cluster_rows = [row for row in read_csv_rows(cluster_csv) if row["filename"] in keep_filenames]
    by_filename: dict[str, list[dict[str, object]]] = {}

    rows_by_filename: dict[str, list[dict[str, str]]] = {}
    for row in cluster_rows:
        rows_by_filename.setdefault(row["filename"], []).append(row)

    for filename, filename_rows in rows_by_filename.items():
        embedding_path = embedding_dir / f"{Path(filename).stem}_embeddings.csv"
        if not embedding_path.exists():
            continue

        labels_needed = {row["label"] for row in filename_rows}
        label_to_bbox = {}
        with embedding_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for emb_row in reader:
                label = emb_row["label"]
                if label not in labels_needed:
                    continue
                label_to_bbox[label] = {
                    "bbox_min_row": as_int(emb_row["bbox_min_row"]),
                    "bbox_min_col": as_int(emb_row["bbox_min_col"]),
                    "bbox_max_row": as_int(emb_row["bbox_max_row"]),
                    "bbox_max_col": as_int(emb_row["bbox_max_col"]),
                }

        for row in filename_rows:
            if row["label"] not in label_to_bbox:
                continue
            bbox = label_to_bbox[row["label"]]
            cluster_id = as_int(row["cluster"])
            by_filename.setdefault(filename, []).append(
                {
                    "label": row["label"],
                    "cluster": cluster_id,
                    "cluster_distance": float(row["cluster_distance"]),
                    **bbox,
                }
            )

    for filename, boxes in by_filename.items():
        boxes.sort(key=lambda box: (box["cluster"], box["cluster_distance"], as_int(box["label"])))
    return by_filename


def parse_shape(value: str) -> tuple[int, int]:
    parts = value.split("x")
    if len(parts) < 2:
        raise ValueError(f"Unexpected shape string: {value}")
    return int(parts[0]), int(parts[1])


def load_image_dims(segmentation_manifest: Path) -> dict[str, dict[str, int]]:
    dims: dict[str, dict[str, int]] = {}
    for row in read_csv_rows(segmentation_manifest):
        source_height, source_width = parse_shape(row["source_image_shape"])
        dims[row["filename"]] = {
            "source_height": source_height,
            "source_width": source_width,
            "crop_height": as_int(row["crop_height"]),
            "crop_width": as_int(row["crop_width"]),
        }
    return dims


def cluster_color(cluster_id: int) -> str:
    return PALETTE[(cluster_id - 1) % len(PALETTE)]


def encode_png_data_uri(path: Path) -> str:
    png_bytes = path.read_bytes()
    encoded = base64.b64encode(png_bytes).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def build_box_svg(
    boxes: list[dict[str, object]],
    source_width: int,
    source_height: int,
    crop_width: int,
    crop_height: int,
) -> str:
    parts = [
        f'<svg class="overlay" viewBox="0 0 {source_width} {source_height}" preserveAspectRatio="none" aria-hidden="true">'
    ]
    parts.append(
        f'<rect class="crop-frame" x="0" y="0" width="{crop_width}" height="{crop_height}"></rect>'
    )
    for box in boxes:
        min_col = int(box["bbox_min_col"])
        min_row = int(box["bbox_min_row"])
        width = max(1, int(box["bbox_max_col"]) - min_col)
        height = max(1, int(box["bbox_max_row"]) - min_row)
        color = cluster_color(int(box["cluster"]))
        parts.append(
            (
                '<rect class="bbox" '
                f'x="{min_col}" y="{min_row}" width="{width}" height="{height}" '
                f'stroke="{color}" data-cluster="{int(box["cluster"])}" data-label="{html.escape(str(box["label"]))}"></rect>'
            )
        )
    parts.append("</svg>")
    return "".join(parts)


def build_legend() -> str:
    items = []
    for cluster_id in range(1, 21):
        color = cluster_color(cluster_id)
        items.append(
            '<div class="legend-item">'
            f'<span class="swatch" style="background:{color}"></span>'
            f"<span>Cluster {cluster_id}</span>"
            "</div>"
        )
    return "".join(items)


def build_card(row: dict[str, object], boxes: list[dict[str, object]]) -> str:
    png_path = Path(row["png_path"])
    img_src = html.escape(encode_png_data_uri(png_path))
    title = f'P{row["passage_position"]}'
    subtitle = f'day {row["time_hours"] / 24:.2f} | {row["passage_id"]}'
    filename = html.escape(str(row["filename"]))
    box_count = len(boxes)
    return (
        '<article class="card">'
        '<div class="image-wrap">'
        f'<img loading="lazy" src="{img_src}" alt="{filename}">'
        f"{build_box_svg(boxes, int(row['source_width']), int(row['source_height']), int(row['crop_width']), int(row['crop_height']))}"
        "</div>"
        '<div class="meta">'
        f"<strong>{html.escape(title)}</strong>"
        f"<div>{html.escape(subtitle)}</div>"
        f'<div class="filename">{filename}</div>'
        f'<div class="counts">{box_count} clustered cells outlined</div>'
        "</div>"
        "</article>"
    )


def build_html(
    lineage_id: str,
    field: str,
    rows: list[dict[str, object]],
    overlays: dict[str, list[dict[str, object]]],
) -> str:
    card_html = []
    total_boxes = 0
    for row in rows:
        boxes = overlays.get(str(row["filename"]), [])
        total_boxes += len(boxes)
        card_html.append(build_card(row, boxes))

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(lineage_id)} stack clusters ({html.escape(field)})</title>
  <style>
    :root {{
      --bg: #f3efe8;
      --panel: #fffdf9;
      --ink: #171512;
      --muted: #655f58;
      --line: #d7cdbf;
      --accent: #8c5a31;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background:
        radial-gradient(circle at top left, rgba(184, 116, 59, 0.10), transparent 28%),
        linear-gradient(180deg, #ece5d8 0%, var(--bg) 100%);
      color: var(--ink);
      font: 15px/1.4 Georgia, "Times New Roman", serif;
    }}
    .page {{
      width: min(1200px, calc(100vw - 32px));
      margin: 20px auto 48px;
    }}
    .intro {{
      position: sticky;
      top: 0;
      z-index: 20;
      display: grid;
      grid-template-columns: 1.2fr 1fr;
      gap: 16px;
      margin-bottom: 18px;
      padding: 16px 18px;
      background: rgba(255, 253, 249, 0.97);
      border: 1px solid var(--line);
      backdrop-filter: blur(6px);
    }}
    .intro h1 {{
      margin: 0 0 6px;
      font-size: 28px;
    }}
    .intro p {{
      margin: 0 0 6px;
      color: var(--muted);
    }}
    .legend {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 8px 10px;
      align-content: start;
    }}
    .legend-item {{
      display: flex;
      align-items: center;
      gap: 8px;
      font-size: 13px;
    }}
    .swatch {{
      width: 16px;
      height: 16px;
      border: 1px solid rgba(0, 0, 0, 0.25);
      flex: 0 0 auto;
    }}
    .controls {{
      display: flex;
      gap: 10px;
      align-items: center;
      margin-top: 10px;
      color: var(--muted);
      font-size: 13px;
    }}
    .stack {{
      display: flex;
      flex-direction: column;
      gap: 14px;
    }}
    .card {{
      background: var(--panel);
      border: 1px solid var(--line);
      overflow: hidden;
      box-shadow: 0 10px 28px rgba(52, 35, 18, 0.08);
    }}
    .image-wrap {{
      position: relative;
      line-height: 0;
      background: #ddd8cf;
    }}
    .image-wrap img {{
      display: block;
      width: 100%;
      height: auto;
    }}
    .overlay {{
      position: absolute;
      inset: 0;
      width: 100%;
      height: 100%;
      pointer-events: none;
    }}
    .bbox {{
      fill: none;
      stroke-width: 1.75;
      vector-effect: non-scaling-stroke;
      opacity: 0.9;
    }}
    .crop-frame {{
      fill: none;
      stroke: rgba(255, 255, 255, 0.6);
      stroke-width: 2.5;
      stroke-dasharray: 8 6;
      vector-effect: non-scaling-stroke;
    }}
    body.hide-boxes .bbox {{
      display: none;
    }}
    body.hide-boxes .crop-frame {{
      display: none;
    }}
    .meta {{
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 8px 12px;
      padding: 10px 12px 12px;
      align-items: start;
    }}
    .meta strong {{
      color: var(--accent);
      font-size: 18px;
    }}
    .filename {{
      color: var(--muted);
      word-break: break-all;
      font-size: 12px;
    }}
    .counts {{
      color: var(--muted);
      font-size: 12px;
      text-align: right;
      white-space: nowrap;
    }}
    @media (max-width: 900px) {{
      .intro {{
        grid-template-columns: 1fr;
      }}
      .legend {{
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }}
      .meta {{
        grid-template-columns: 1fr;
      }}
      .counts {{
        text-align: left;
      }}
    }}
  </style>
</head>
<body>
  <main class="page">
    <header class="intro">
      <div>
        <h1>{html.escape(lineage_id)}</h1>
        <p>{len(rows)} cached PNGs in field <strong>{html.escape(field)}</strong>.</p>
        <p>{total_boxes} clustered cells outlined using pooled k-means cluster assignments.</p>
        <div class="controls">
          <label><input id="toggle-boxes" type="checkbox" checked> Show bounding boxes</label>
        </div>
      </div>
      <div>
        <div class="legend">{build_legend()}</div>
      </div>
    </header>
    <section class="stack">
      {''.join(card_html)}
    </section>
  </main>
  <script>
    const toggle = document.getElementById("toggle-boxes");
    toggle.addEventListener("change", () => {{
      document.body.classList.toggle("hide-boxes", !toggle.checked);
    }});
  </script>
</body>
</html>
"""


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    lineage_name = lineage_dir_name(args.lineage_id)

    image_summary_path = resolve_path(
        project_root,
        args.image_summary or f"data/lineage_area/{lineage_name}/image_summary.csv",
    )
    cluster_csv = resolve_path(
        project_root,
        args.cluster_csv or f"data/lineage_embedding_clusters/{lineage_name}/pooled_all_quantiles/best_clustering.csv",
    )
    segmentation_manifest = resolve_path(project_root, args.segmentation_manifest)
    embedding_dir = resolve_path(project_root, args.embedding_dir)
    png_dir = resolve_path(project_root, args.png_dir or f"data/lineage_png/{lineage_name}")
    out_html = resolve_path(project_root, args.out_html)

    if image_summary_path is None or not image_summary_path.exists():
      raise FileNotFoundError(f"Missing image summary: {image_summary_path}")
    if cluster_csv is None or not cluster_csv.exists():
      raise FileNotFoundError(f"Missing cluster CSV: {cluster_csv}")
    if segmentation_manifest is None or not segmentation_manifest.exists():
      raise FileNotFoundError(f"Missing segmentation manifest: {segmentation_manifest}")
    if embedding_dir is None or not embedding_dir.exists():
      raise FileNotFoundError(f"Missing embedding directory: {embedding_dir}")
    if png_dir is None or not png_dir.exists():
      raise FileNotFoundError(f"Missing PNG directory: {png_dir}")
    if out_html is None:
      raise ValueError("Output HTML path is required")

    image_dims = load_image_dims(segmentation_manifest)
    rows = load_stack_rows(image_summary_path, png_dir, args.field, image_dims)
    keep_filenames = {str(row["filename"]) for row in rows}
    overlays = load_cluster_boxes(cluster_csv, embedding_dir, keep_filenames)

    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(
        build_html(args.lineage_id, args.field, rows, overlays),
        encoding="utf-8",
    )
    print(f"Wrote {out_html}")


if __name__ == "__main__":
    main()
