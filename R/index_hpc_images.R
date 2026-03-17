parse_args <- function(args) {
  parsed <- list(
    passaging_csv = "core_data/passaging.csv",
    media_csv = "core_data/media.csv",
    out_csv = "data/image_manifest.csv",
    out_id_csv = "data/image_id_summary.csv"
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

  required <- c("hpc_dir", "passaging_csv", "media_csv", "out_csv", "out_id_csv")
  missing <- required[!nzchar(vapply(required, function(x) {
    val <- parsed[[x]]
    if (is.null(val)) "" else as.character(val)
  }, character(1)))]
  if (length(missing) > 0) {
    stop(sprintf("Missing required args: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  parsed
}

ensure_parent_dir <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
}

parse_filename <- function(filename) {
  m <- regexec("^(.*)_10x_ph_(bl|br|tl|tr)(?:_(overlay))?\\.([A-Za-z0-9]+)$", filename)
  parts <- regmatches(filename, m)[[1]]

  if (length(parts) == 0) {
    return(data.frame(
      filename = filename,
      parse_ok = FALSE,
      id = NA_character_,
      field = NA_character_,
      overlay_tag = NA_character_,
      ext = tolower(sub(".*\\.([^.]+)$", "\\1", filename)),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    filename = filename,
    parse_ok = TRUE,
    id = parts[[2]],
    field = parts[[3]],
    overlay_tag = if (length(parts) >= 4 && nzchar(parts[[4]])) parts[[4]] else "",
    ext = tolower(parts[[5]]),
    stringsAsFactors = FALSE
  )
}

artifact_type_from_parts <- function(ext, overlay_tag) {
  if (!is.na(overlay_tag) && nzchar(overlay_tag)) {
    return("overlay")
  }
  if (tolower(ext) %in% c("tif", "tiff")) {
    return("raw_image")
  }
  if (tolower(ext) == "png") {
    return("png_sidecar")
  }
  "other"
}

derive_condition_label <- function(stressor, concentration, unit) {
  stressor <- ifelse(is.na(stressor), "", trimws(stressor))
  concentration <- ifelse(is.na(concentration), "", trimws(as.character(concentration)))
  unit <- ifelse(is.na(unit), "", trimws(unit))

  out <- rep("untreated", length(stressor))
  has_stressor <- nzchar(stressor)
  out[has_stressor] <- trimws(paste(stressor[has_stressor], concentration[has_stressor], unit[has_stressor]))
  out[has_stressor] <- gsub("\\s+", " ", out[has_stressor])
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

files <- list.files(args$hpc_dir, recursive = FALSE, full.names = TRUE)
files <- files[file.info(files)$isdir %in% FALSE]

manifest <- do.call(
  rbind,
  lapply(basename(files), parse_filename)
)

manifest$filepath <- normalizePath(files, winslash = "/", mustWork = FALSE)
manifest$artifact_type <- mapply(artifact_type_from_parts, manifest$ext, manifest$overlay_tag, USE.NAMES = FALSE)
manifest$tile_group <- ifelse(manifest$parse_ok, paste0("10x_ph_", manifest$field), NA_character_)

passaging <- read.csv(args$passaging_csv, stringsAsFactors = FALSE, check.names = FALSE)
media <- read.csv(args$media_csv, stringsAsFactors = FALSE, check.names = FALSE)

manifest <- merge(manifest, passaging, by = "id", all.x = TRUE, sort = FALSE)
manifest <- merge(
  manifest,
  media,
  by.x = "media",
  by.y = "id",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_media")
)

manifest$matched_passaging <- !is.na(manifest$cellLine)
manifest$condition_label <- derive_condition_label(
  manifest$Stressor,
  manifest$Stressor_concentration,
  manifest$Stressor_unit
)

id_summary <- unique(manifest[c(
  "id", "cellLine", "event", "passage", "date", "address", "comment", "media",
  "growthType", "matched_passaging", "condition_label"
)])

raw_counts <- aggregate(
  list(raw_file_count = manifest$artifact_type == "raw_image"),
  by = list(id = manifest$id),
  FUN = sum
)
overlay_counts <- aggregate(
  list(overlay_file_count = manifest$artifact_type == "overlay"),
  by = list(id = manifest$id),
  FUN = sum
)

id_summary <- merge(id_summary, raw_counts, by = "id", all.x = TRUE, sort = FALSE)
id_summary <- merge(id_summary, overlay_counts, by = "id", all.x = TRUE, sort = FALSE)
id_summary$raw_file_count[is.na(id_summary$raw_file_count)] <- 0L
id_summary$overlay_file_count[is.na(id_summary$overlay_file_count)] <- 0L
id_summary$has_complete_raw_tiles <- id_summary$raw_file_count == 4L

ensure_parent_dir(args$out_csv)
ensure_parent_dir(args$out_id_csv)

write.csv(manifest, args$out_csv, row.names = FALSE, na = "")
write.csv(id_summary, args$out_id_csv, row.names = FALSE, na = "")

cat(sprintf("Indexed %s files from %s\n", nrow(manifest), args$hpc_dir))
cat(sprintf("Parsed filenames: %s\n", sum(manifest$parse_ok)))
cat(sprintf("Matched passaging IDs: %s\n", sum(id_summary$matched_passaging, na.rm = TRUE)))
cat(sprintf("Raw images: %s\n", sum(manifest$artifact_type == "raw_image")))
cat(sprintf("Overlays: %s\n", sum(manifest$artifact_type == "overlay")))
cat(sprintf("Wrote manifest: %s\n", args$out_csv))
cat(sprintf("Wrote ID summary: %s\n", args$out_id_csv))
