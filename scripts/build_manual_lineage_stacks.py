#!/usr/bin/env python3
import argparse
import base64
import csv
import html
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build one standalone HTML stack per manually curated lineage using cached PNGs."
    )
    parser.add_argument(
        "--project-root",
        default=".",
        help="Project root containing data/lineage_area and data/lineage_png",
    )
    parser.add_argument(
        "--field",
        default="tl",
        help="Image field to include in the stack (for example tl, tr, bl, br)",
    )
    parser.add_argument(
        "--out-dir",
        default="data/lineage_png/stacks",
        help="Output directory for HTML files, relative to project root unless absolute",
    )
    return parser.parse_args()


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def lineage_sort_key(name: str) -> tuple[int, str]:
    suffix = name.removeprefix("lineage_manual_").removeprefix("manual_")
    try:
        return (int(suffix), name)
    except ValueError:
        return (10**9, name)


def load_lineage_rows(project_root: Path, field: str) -> dict[str, list[dict[str, object]]]:
    area_root = project_root / "data" / "lineage_area"
    png_root = project_root / "data" / "lineage_png"
    out: dict[str, list[dict[str, object]]] = {}

    for lineage_dir in sorted(area_root.glob("lineage_manual_*"), key=lambda p: lineage_sort_key(p.name)):
        summary_path = lineage_dir / "image_summary.csv"
        if not summary_path.exists():
            continue

        lineage_id = lineage_dir.name.removeprefix("lineage_")
        rows: list[dict[str, object]] = []
        with summary_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                if row.get("field") != field:
                    continue

                filename = row["filename"]
                png_name = Path(filename).with_suffix(".png").name
                png_path = png_root / lineage_dir.name / png_name
                if not png_path.exists():
                    continue

                rows.append(
                    {
                        "lineage_id": lineage_id,
                        "passage_position": int(row["passage_position"]),
                        "passage_id": row["passage_id"],
                        "png_name": png_name,
                        "time_hours": as_float(row["mean_time_since_passage_start_hours"]),
                        "png_path": png_path,
                    }
                )

        rows.sort(key=lambda row: (row["passage_position"], row["time_hours"], row["png_name"]))
        out[lineage_id] = rows

    return out


def encode_png_data_uri(path: Path) -> str:
    png_bytes = path.read_bytes()
    encoded = base64.b64encode(png_bytes).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def build_card(row: dict[str, object]) -> str:
    time_hours = float(row["time_hours"])
    time_days = time_hours / 24 if time_hours == time_hours else float("nan")
    img_src = html.escape(encode_png_data_uri(Path(row["png_path"])))
    title = f'{row["lineage_id"]} | P{row["passage_position"]}'
    subtitle = f'day {time_days:.2f} | {html.escape(str(row["passage_id"]))}'
    png_label = html.escape(str(row["png_name"]))
    return (
        '<article class="card">'
        f'<img loading="lazy" src="{img_src}" alt="{html.escape(title)}">'
        '<div class="meta">'
        f"<strong>{html.escape(title)}</strong>"
        f"<div>{subtitle}</div>"
        f'<div class="filename">{png_label}</div>'
        "</div>"
        "</article>"
    )


def build_html(lineage_id: str, rows: list[dict[str, object]], field: str) -> str:
    cards = "".join(build_card(row) for row in rows)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(lineage_id)} stack ({html.escape(field)})</title>
  <style>
    :root {{
      --bg: #f4f1ea;
      --panel: #fffdf8;
      --ink: #1f1d1a;
      --muted: #6b655d;
      --line: #d8d0c4;
      --accent: #7b4b2a;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: linear-gradient(180deg, #ede7dc 0%, var(--bg) 100%);
      color: var(--ink);
      font: 15px/1.4 Georgia, "Times New Roman", serif;
    }}
    .page {{
      width: min(920px, calc(100vw - 32px));
      margin: 20px auto 48px;
    }}
    .intro {{
      position: sticky;
      top: 0;
      z-index: 10;
      margin-bottom: 18px;
      padding: 16px 18px;
      background: rgba(255, 253, 248, 0.96);
      border: 1px solid var(--line);
    }}
    .intro h1 {{
      margin: 0 0 4px;
      font-size: 28px;
    }}
    .intro p {{
      margin: 0;
      color: var(--muted);
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
      box-shadow: 0 8px 20px rgba(54, 37, 20, 0.08);
    }}
    .card img {{
      display: block;
      width: 100%;
      height: auto;
      background: #ddd8cf;
    }}
    .meta {{
      padding: 10px 12px 12px;
    }}
    .meta strong {{
      color: var(--accent);
    }}
    .filename {{
      color: var(--muted);
      word-break: break-all;
      font-size: 12px;
      margin-top: 2px;
    }}
  </style>
</head>
<body>
  <main class="page">
    <header class="intro">
      <h1>{html.escape(lineage_id)}</h1>
      <p>{len(rows)} cached PNGs, ordered by passage then time since seeding. Field: <strong>{html.escape(field)}</strong>.</p>
    </header>
    <section class="stack">{cards}</section>
  </main>
</body>
</html>
"""


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = project_root / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    lineage_rows = load_lineage_rows(project_root, args.field)
    for lineage_id, rows in sorted(lineage_rows.items(), key=lambda item: lineage_sort_key(item[0])):
        html_path = out_dir / f"{lineage_id}_{args.field}_stack.html"
        html_path.write_text(build_html(lineage_id, rows, args.field), encoding="utf-8")
        print(f"Wrote {html_path}")


if __name__ == "__main__":
    main()
