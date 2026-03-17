parse_args <- function(args) {
  parsed <- list(
    out_image_csv = "manifests/dev_subset_images.csv",
    out_id_csv = "manifests/dev_subset_ids.csv",
    copy_script = "scripts/copy_dev_subset.sh",
    n_groups = "4",
    ids_per_group = "4"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for %s", key), call. = FALSE)
    }
    parsed[[substring(key, 3)]] <- args[[i + 1]]
    i <- i + 2
  }

  required <- c("manifest_csv", "out_image_csv", "out_id_csv", "copy_script")
  missing <- required[!nzchar(vapply(required, function(x) {
    val <- parsed[[x]]
    if (is.null(val)) "" else as.character(val)
  }, character(1)))]
  if (length(missing) > 0) {
    stop(sprintf("Missing required args: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  parsed$n_groups <- as.integer(parsed$n_groups)
  parsed$ids_per_group <- as.integer(parsed$ids_per_group)
  parsed
}

ensure_parent_dir <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
}

coalesce_chr <- function(x, y) {
  ifelse(is.na(x) | !nzchar(trimws(x)), y, x)
}

pick_spread_indices <- function(n, k) {
  if (n <= 0 || k <= 0) {
    return(integer(0))
  }
  if (n <= k) {
    return(seq_len(n))
  }

  idx <- unique(round(seq(1, n, length.out = k)))
  idx <- pmax(1, pmin(n, idx))
  while (length(idx) < k) {
    remaining <- setdiff(seq_len(n), idx)
    idx <- sort(c(idx, remaining[[1]]))
  }
  idx
}

build_copy_script <- function(filepaths, script_path) {
  lines <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "",
    'SRC_ROOT="${1:-/share/lab_crd/lab_crd/CLONEID/data/LTEEs}"',
    'DEST_ROOT="${2:-./dev_data/raw}"',
    'mkdir -p "${DEST_ROOT}"',
    ""
  )

  relpaths <- basename(filepaths)
  copy_lines <- sprintf('cp -n "${SRC_ROOT}/%s" "${DEST_ROOT}/"', relpaths)
  lines <- c(lines, copy_lines)
  writeLines(lines, script_path, useBytes = TRUE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
manifest <- read.csv(args$manifest_csv, stringsAsFactors = FALSE, check.names = FALSE)

raw <- manifest[
  manifest$artifact_type == "raw_image" &
    manifest$parse_ok &
    manifest$matched_passaging,
]

if (nrow(raw) == 0) {
  stop("No matched raw images found in manifest.", call. = FALSE)
}

raw$id <- as.character(raw$id)
raw$cellLine <- coalesce_chr(raw$cellLine, "unknown_cell_line")
raw$condition_label <- coalesce_chr(raw$condition_label, "untreated")
raw$event <- coalesce_chr(raw$event, "unknown_event")
raw$date <- coalesce_chr(raw$date, "")
raw$passage_num <- suppressWarnings(as.numeric(raw$passage))
raw$field <- coalesce_chr(raw$field, "unknown_field")

id_level <- unique(raw[c(
  "id", "cellLine", "condition_label", "event", "passage", "passage_num",
  "date", "growthType", "media", "address"
)])

field_counts <- aggregate(list(field_count = raw$field), by = list(id = raw$id), FUN = function(x) length(unique(x)))
id_level <- merge(id_level, field_counts, by = "id", all.x = TRUE, sort = FALSE)
id_level$field_count[is.na(id_level$field_count)] <- 0L
id_level <- id_level[id_level$field_count >= 4, ]

if (nrow(id_level) == 0) {
  stop("No IDs with four raw-image tiles were found.", call. = FALSE)
}

id_level$group_key <- paste(id_level$cellLine, id_level$condition_label, sep = " | ")
group_sizes <- aggregate(list(id_count = id_level$id), by = list(group_key = id_level$group_key), FUN = length)
group_sizes <- group_sizes[order(-group_sizes$id_count, group_sizes$group_key), ]
selected_groups <- head(group_sizes$group_key, args$n_groups)

selected_ids <- character(0)
for (group in selected_groups) {
  group_rows <- id_level[id_level$group_key == group, ]
  group_rows <- group_rows[order(
    is.na(group_rows$passage_num),
    group_rows$passage_num,
    group_rows$date,
    group_rows$id
  ), ]
  picks <- pick_spread_indices(nrow(group_rows), args$ids_per_group)
  selected_ids <- c(selected_ids, group_rows$id[picks])
}

selected_ids <- unique(selected_ids)
selected_id_level <- id_level[id_level$id %in% selected_ids, ]
selected_id_level <- selected_id_level[order(selected_id_level$group_key, selected_id_level$passage_num, selected_id_level$date, selected_id_level$id), ]

selected_images <- raw[raw$id %in% selected_ids, ]
selected_images <- merge(
  selected_images,
  selected_id_level[c("id", "group_key")],
  by = "id",
  all.x = TRUE,
  sort = FALSE
)
selected_images <- selected_images[order(selected_images$group_key, selected_images$id, selected_images$field), ]

ensure_parent_dir(args$out_image_csv)
ensure_parent_dir(args$out_id_csv)
ensure_parent_dir(args$copy_script)

write.csv(selected_images, args$out_image_csv, row.names = FALSE, na = "")
write.csv(selected_id_level, args$out_id_csv, row.names = FALSE, na = "")
build_copy_script(selected_images$filepath, args$copy_script)

cat(sprintf("Selected groups: %s\n", paste(selected_groups, collapse = "; ")))
cat(sprintf("Selected IDs: %s\n", length(unique(selected_images$id))))
cat(sprintf("Selected raw images: %s\n", nrow(selected_images)))
cat(sprintf("Wrote image subset: %s\n", args$out_image_csv))
cat(sprintf("Wrote ID subset: %s\n", args$out_id_csv))
cat(sprintf("Wrote copy script: %s\n", args$copy_script))
