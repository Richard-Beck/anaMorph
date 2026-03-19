parse_args <- function(args) {
  parsed <- list(
    segmentation_manifest = "manifests/dev_subset_segmentation_images.csv",
    metadata_manifest = "manifests/dev_subset_images.csv",
    passaging_csv = "core_data/passaging.csv",
    output_dir = "data/dev_unsupervised_analysis",
    embedding_pcs = "20"
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

  parsed$embedding_pcs <- as.integer(parsed$embedding_pcs)
  parsed
}

resolve_project_root <- function() {
  candidates <- c(
    normalizePath(".", winslash = "/", mustWork = FALSE),
    normalizePath("..", winslash = "/", mustWork = FALSE),
    normalizePath("../..", winslash = "/", mustWork = FALSE)
  )

  for (candidate in unique(candidates)) {
    if (file.exists(file.path(candidate, "README.md")) &&
        dir.exists(file.path(candidate, "manifests")) &&
        dir.exists(file.path(candidate, "dev_data"))) {
      return(candidate)
    }
  }

  stop("Could not locate project root from current working directory.", call. = FALSE)
}

project_path <- function(project_root, path) {
  is_absolute <- grepl("^[A-Za-z]:[/\\\\]|^/", path)
  out <- path
  out[!is_absolute] <- file.path(project_root, path[!is_absolute])
  out
}

normalize_relpath <- function(path) {
  gsub("\\\\", "/", path)
}

read_mask_tiff <- function(path) {
  tiff::readTIFF(path, as.is = TRUE)
}

read_grayscale_image <- function(path) {
  img <- tiff::readTIFF(path, native = FALSE)
  dims <- dim(img)
  if (is.null(dims) || length(dims) == 2) {
    return(img)
  }
  if (length(dims) == 3) {
    channels <- seq_len(min(3, dims[[3]]))
    return(apply(img[, , channels, drop = FALSE], c(1, 2), mean))
  }
  stop(sprintf("Unsupported image dimensionality for %s", path), call. = FALSE)
}

derive_mask_features <- function(seg_row, embedding_rows) {
  mask <- read_mask_tiff(seg_row$mask_relpath)
  image <- read_grayscale_image(seg_row$raw_relpath)
  mask_height <- nrow(mask)
  mask_width <- ncol(mask)

  purrr::pmap_dfr(
    list(
      embedding_rows$label,
      embedding_rows$bbox_min_row,
      embedding_rows$bbox_min_col,
      embedding_rows$bbox_max_row,
      embedding_rows$bbox_max_col,
      embedding_rows$area_px
    ),
    function(label, min_row, min_col, max_row, max_col, area_px) {
      row_idx <- seq.int(min_row + 1L, max_row)
      col_idx <- seq.int(min_col + 1L, max_col)
      object_mask <- mask[row_idx, col_idx, drop = FALSE] == label
      object_image <- image[row_idx, col_idx, drop = FALSE]
      object_image[!object_mask] <- 0

      bbox_height <- max_row - min_row
      bbox_width <- max_col - min_col
      bbox_area <- bbox_height * bbox_width
      object_pixels <- which(object_mask, arr.ind = TRUE)
      centroid_row <- mean(object_pixels[, 1]) + min_row
      centroid_col <- mean(object_pixels[, 2]) + min_col

      tibble::tibble(
        label = label,
        bbox_height = bbox_height,
        bbox_width = bbox_width,
        bbox_area = bbox_area,
        bbox_aspect = bbox_height / pmax(bbox_width, 1),
        bbox_fill_fraction = area_px / pmax(bbox_area, 1),
        edge_touching = min_row == 0L | min_col == 0L | max_row >= mask_height | max_col >= mask_width,
        centroid_row = centroid_row,
        centroid_col = centroid_col,
        centroid_row_scaled = centroid_row / mask_height,
        centroid_col_scaled = centroid_col / mask_width,
        intensity_mean = mean(object_image[object_mask]),
        intensity_sd = stats::sd(object_image[object_mask]),
        intensity_q90 = stats::quantile(object_image[object_mask], probs = 0.9, names = FALSE),
        image_height = mask_height,
        image_width = mask_width
      )
    }
  )
}

main <- function() {
  required_packages <- c("dplyr", "purrr", "readr", "tibble", "tiff")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      sprintf("Install required packages before running: %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }

  args <- parse_args(commandArgs(trailingOnly = TRUE))
  project_root <- resolve_project_root()

  segmentation_manifest <- project_path(project_root, args$segmentation_manifest)
  metadata_manifest <- project_path(project_root, args$metadata_manifest)
  passaging_csv <- project_path(project_root, args$passaging_csv)
  output_dir <- project_path(project_root, args$output_dir)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  seg_manifest <- readr::read_csv(segmentation_manifest, show_col_types = FALSE) |>
    dplyr::mutate(
      raw_relpath = project_path(project_root, normalize_relpath(raw_relpath)),
      mask_relpath = project_path(project_root, normalize_relpath(mask_relpath)),
      embedding_relpath = project_path(project_root, normalize_relpath(embedding_relpath))
    )

  passaging_metadata <- readr::read_csv(passaging_csv, show_col_types = FALSE) |>
    dplyr::select(id, event) |>
    dplyr::distinct()

  seg_manifest <- seg_manifest |>
    dplyr::left_join(passaging_metadata, by = "id") |>
    dplyr::filter(is.na(event) | event != "seeding")

  dev_metadata <- readr::read_csv(metadata_manifest, show_col_types = FALSE) |>
    dplyr::select(
      filename, id, field, cellLine, event, passage,
      passage_num, date, growthType, media, group_key
    ) |>
    dplyr::distinct()

  embedding_tbl <- purrr::map_dfr(seq_len(nrow(seg_manifest)), function(i) {
    seg_row <- seg_manifest[i, ]
    readr::read_csv(seg_row$embedding_relpath, show_col_types = FALSE) |>
      dplyr::mutate(
        id = seg_row$id,
        field = seg_row$field,
        filename = seg_row$filename,
        raw_relpath = seg_row$raw_relpath,
        mask_relpath = seg_row$mask_relpath,
        run_id = seg_row$run_id,
        timestamp_utc = seg_row$timestamp_utc
      )
  })

  mask_feature_tbl <- purrr::map_dfr(seq_len(nrow(seg_manifest)), function(i) {
    seg_row <- seg_manifest[i, ]
    embedding_rows <- embedding_tbl |>
      dplyr::filter(filename == seg_row$filename)
    derive_mask_features(seg_row, embedding_rows) |>
      dplyr::mutate(filename = seg_row$filename)
  })

  object_tbl <- embedding_tbl |>
    dplyr::left_join(mask_feature_tbl, by = c("filename", "label")) |>
    dplyr::left_join(dev_metadata, by = c("filename", "id", "field")) |>
    dplyr::mutate(
      cellLine = dplyr::coalesce(cellLine, "unknown"),
      edge_touching = as.logical(edge_touching)
    )

  embedding_features <- names(object_tbl)[grepl("^embedding_", names(object_tbl))]
  embedding_mat <- object_tbl[, embedding_features, drop = FALSE] |> as.matrix()
  keep_rows <- stats::complete.cases(embedding_mat)
  embedding_scaled <- scale(embedding_mat[keep_rows, , drop = FALSE])
  embedding_pca <- stats::prcomp(embedding_scaled, center = FALSE, scale. = FALSE)
  n_pcs <- min(args$embedding_pcs, ncol(embedding_pca$x))

  embedding_pca_tbl <- object_tbl[keep_rows, c("filename", "label", "id", "field", "cellLine"), drop = FALSE] |>
    dplyr::bind_cols(
      tibble::as_tibble(embedding_pca$x[, seq_len(n_pcs), drop = FALSE], .name_repair = "minimal") |>
        stats::setNames(sprintf("embedding_PC%02d", seq_len(n_pcs)))
    )

  saveRDS(seg_manifest, file.path(output_dir, "seg_manifest_filtered.rds"))
  saveRDS(embedding_tbl, file.path(output_dir, "embedding_table.rds"))
  saveRDS(mask_feature_tbl, file.path(output_dir, "mask_feature_table.rds"))
  saveRDS(object_tbl, file.path(output_dir, "dev_object_table.rds"))
  saveRDS(embedding_pca, file.path(output_dir, "embedding_pca_model.rds"))
  saveRDS(embedding_pca_tbl, file.path(output_dir, "embedding_pca_table.rds"))

  readr::write_csv(object_tbl, file.path(output_dir, "dev_object_table.csv"))
  readr::write_csv(embedding_pca_tbl, file.path(output_dir, "embedding_pca_table.csv"))

  cat(sprintf("Filtered images: %d\n", nrow(seg_manifest)))
  cat(sprintf("Objects retained: %d\n", nrow(object_tbl)))
  cat(sprintf("Saved preprocessing outputs to %s\n", normalize_relpath(output_dir)))
}

main()
