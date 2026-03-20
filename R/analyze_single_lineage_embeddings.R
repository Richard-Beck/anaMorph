parse_args <- function(args) {
  parsed <- list(
    lineage_rds = "core_data/lineages.Rds",
    segmentation_manifest = "data/all_images/manifests/segmentation_images.csv",
    embedding_dir = "data/all_images/embeddings",
    out_dir = "data/lineage_embeddings",
    n_pcs = "20",
    top_corr_pcs = "3",
    top_export_pcs = "5",
    crops_per_side = "10",
    n_cores = "1"
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

  if (is.null(parsed$lineage_id) || !nzchar(parsed$lineage_id)) {
    stop("Missing required arg: --lineage_id", call. = FALSE)
  }

  parsed$n_pcs <- as.integer(parsed$n_pcs)
  parsed$top_corr_pcs <- as.integer(parsed$top_corr_pcs)
  parsed$top_export_pcs <- as.integer(parsed$top_export_pcs)
  parsed$crops_per_side <- as.integer(parsed$crops_per_side)
  parsed$n_cores <- as.integer(parsed$n_cores)
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
        dir.exists(file.path(candidate, "core_data"))) {
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

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

warning_state <- new.env(parent = emptyenv())
warning_state$messages <- character()
warning_state$log_path <- NULL

reset_warning_state <- function(log_path = NULL) {
  warning_state$messages <- character()
  warning_state$log_path <- log_path
  if (!is.null(log_path) && file.exists(log_path)) {
    file.remove(log_path)
  }
}

emit_warning <- function(message_text) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] WARNING: %s", timestamp, message_text)
  warning_state$messages <- c(warning_state$messages, line)
  cat(line, "\n", sep = "")
  if (!is.null(warning_state$log_path)) {
    write(line, file = warning_state$log_path, append = TRUE)
  }
}

extract_lineage <- function(lineages, lineage_id) {
  if (lineage_id %in% names(lineages)) {
    return(lineages[[lineage_id]])
  }

  hit <- vapply(
    lineages,
    function(x) identical(as.character(x$lineage_id), as.character(lineage_id)),
    logical(1)
  )
  if (!any(hit)) {
    stop(sprintf("Lineage id not found: %s", lineage_id), call. = FALSE)
  }

  lineages[[which(hit)[[1]]]]
}

build_lineage_image_meta <- function(lineage_obj) {
  passage_positions <- stats::setNames(seq_along(lineage_obj$passage_ids), lineage_obj$passage_ids)
  meta_list <- lapply(names(lineage_obj$passages), function(passage_id) {
    current_passage_id <- as.character(passage_id)
    current_passage_number <- as.integer(passage_positions[[current_passage_id]])
    passage_dt <- data.table::as.data.table(lineage_obj$passages[[current_passage_id]])
    if (!nrow(passage_dt)) {
      return(NULL)
    }

    out <- data.table::copy(passage_dt)
    keep_cols <- intersect(c("id", "filename", "filepath", "field", "date"), names(out))
    keep_dt <- out[, ..keep_cols]
    keep_dt[, image_id := if ("filename" %in% names(keep_dt)) {
      ifelse(!is.na(filename) & nzchar(filename), as.character(filename), as.character(id))
    } else {
      as.character(id)
    }]
    keep_dt[, passage_id := current_passage_id]
    keep_dt[, passage_number := current_passage_number]
    keep_dt
  })

  meta_dt <- data.table::rbindlist(meta_list, use.names = TRUE, fill = TRUE)
  if (!nrow(meta_dt)) {
    stop(sprintf("Lineage %s does not contain any images.", lineage_obj$lineage_id), call. = FALSE)
  }

  if (!"id" %in% names(meta_dt)) {
    stop("Lineage image metadata is missing `id`.", call. = FALSE)
  }
  meta_dt[, id := as.character(id)]
  if ("image_id" %in% names(meta_dt)) {
    meta_dt[, image_id := as.character(image_id)]
  }
  if ("filename" %in% names(meta_dt)) {
    meta_dt[, filename := as.character(filename)]
  }
  unique(meta_dt, by = intersect(c("image_id", "filename", "passage_id"), names(meta_dt)))
}

candidate_embedding_paths <- function(project_root, embedding_dir, image_id, filename = NA_character_) {
  stem_candidates <- unique(stats::na.omit(c(
    as.character(image_id),
    tools::file_path_sans_ext(basename(as.character(filename)))
  )))
  file.path(embedding_dir, paste0(stem_candidates, "_embeddings.csv"))
}

match_segmentation_row <- function(seg_dt, row) {
  seg_row <- NULL

  if ("filename" %in% names(seg_dt) && "filename" %in% names(row)) {
    seg_row <- seg_dt[filename == row$filename[[1]]]
    if (nrow(seg_row)) {
      return(seg_row[1])
    }
  }

  if ("source_filepath" %in% names(seg_dt) && "filepath" %in% names(row)) {
    seg_row <- seg_dt[source_filepath == row$filepath[[1]]]
    if (nrow(seg_row)) {
      return(seg_row[1])
    }
  }

  if ("source_basename" %in% names(seg_dt) && "filepath" %in% names(row)) {
    seg_row <- seg_dt[source_basename == basename(row$filepath[[1]])]
    if (nrow(seg_row)) {
      return(seg_row[1])
    }
  }

  if ("id" %in% names(seg_dt)) {
    seg_row <- seg_dt[id == row$id[[1]]]
    if (nrow(seg_row)) {
      return(seg_row[1])
    }
  }

  NULL
}

read_grayscale_image <- function(path) {
  img <- tiff::readTIFF(path, native = FALSE, all = FALSE)
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

read_mask_tiff <- function(path) {
  tiff::readTIFF(path, as.is = TRUE)
}

rescale_crop <- function(x) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  hi <- stats::quantile(x, probs = 0.995, na.rm = TRUE, names = FALSE)
  lo <- min(x, na.rm = TRUE)
  if (!is.finite(hi) || hi <= lo) {
    hi <- max(x, na.rm = TRUE)
  }
  if (!is.finite(hi) || hi <= lo) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x)))
  }
  out <- (x - lo) / (hi - lo)
  out[out < 0] <- 0
  out[out > 1] <- 1
  out
}

resize_nn <- function(mat, out_height = 96, out_width = 96) {
  y_idx <- round(seq(1, nrow(mat), length.out = out_height))
  x_idx <- round(seq(1, ncol(mat), length.out = out_width))
  mat[y_idx, x_idx, drop = FALSE]
}

pad_crop <- function(mat, target_height, target_width) {
  out <- matrix(0, nrow = target_height, ncol = target_width)
  row_offset <- floor((target_height - nrow(mat)) / 2) + 1L
  col_offset <- floor((target_width - ncol(mat)) / 2) + 1L
  out[
    row_offset:(row_offset + nrow(mat) - 1L),
    col_offset:(col_offset + ncol(mat) - 1L)
  ] <- mat
  out
}

assemble_crop_panel <- function(crop_mats, ncol_panel = NULL) {
  if (!length(crop_mats)) {
    return(NULL)
  }
  heights <- vapply(crop_mats, nrow, integer(1))
  widths <- vapply(crop_mats, ncol, integer(1))
  target_height <- max(heights)
  target_width <- max(widths)
  padded <- lapply(crop_mats, pad_crop, target_height = target_height, target_width = target_width)

  n_tiles <- length(padded)
  if (is.null(ncol_panel)) {
    ncol_panel <- ceiling(sqrt(n_tiles))
  }
  nrow_panel <- ceiling(n_tiles / ncol_panel)
  panel <- matrix(0, nrow = nrow_panel * target_height, ncol = ncol_panel * target_width)

  for (i in seq_along(padded)) {
    row_block <- ((i - 1L) %/% ncol_panel)
    col_block <- ((i - 1L) %% ncol_panel)
    row_start <- row_block * target_height + 1L
    col_start <- col_block * target_width + 1L
    panel[
      row_start:(row_start + target_height - 1L),
      col_start:(col_start + target_width - 1L)
    ] <- padded[[i]]
  }

  panel
}

write_crop_tiff <- function(mat, out_path) {
  ensure_dir(dirname(out_path))
  tiff::writeTIFF(rescale_crop(mat), out_path, bits.per.sample = 8L, compression = "LZW")
}

extract_object_crop <- function(image, mask, label, min_row, min_col, max_row, max_col) {
  row_idx <- seq.int(as.integer(min_row) + 1L, as.integer(max_row))
  col_idx <- seq.int(as.integer(min_col) + 1L, as.integer(max_col))
  crop <- image[row_idx, col_idx, drop = FALSE]
  crop_mask <- mask[row_idx, col_idx, drop = FALSE] == as.integer(label)
  crop[!crop_mask] <- 0
  crop
}

assign_count_quantile <- function(cell_count, cutpoints) {
  if (!is.finite(cell_count)) {
    return(NA_character_)
  }
  if (cell_count <= cutpoints[[1]]) {
    return("Q1_low")
  }
  if (cell_count <= cutpoints[[2]]) {
    return("Q2_mid")
  }
  "Q3_high"
}

build_image_count_plot <- function(image_counts, cutpoints, out_path, lineage_id) {
  plt <- ggplot2::ggplot(
    image_counts[!is.na(cell_count)],
    ggplot2::aes(x = cell_count, fill = quantile)
  ) +
    ggplot2::geom_histogram(bins = min(30L, max(10L, nrow(image_counts))), color = "white", alpha = 0.9) +
    ggplot2::geom_vline(
      xintercept = cutpoints,
      linetype = "dashed",
      linewidth = 0.7,
      color = c("#1b9e77", "#d95f02")
    ) +
    ggplot2::scale_fill_manual(
      values = c(Q1_low = "#440154", Q2_mid = "#21908C", Q3_high = "#FDE725"),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = sprintf("Per-image cell count distribution: %s", lineage_id),
      x = "Cell count per image",
      y = "Images",
      fill = "Quantile"
    ) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(out_path, plt, width = 8, height = 5, dpi = 150)
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) {
    return(NA_real_)
  }
  stats::cor(x[ok], y[ok], method = "pearson")
}

build_score_table <- function(object_dt, pca_obj, n_pcs_keep, n_pcs_pad) {
  score_cols <- sprintf("PC_%d", seq_len(n_pcs_keep))
  score_dt <- data.table::copy(object_dt)
  score_values <- data.table::as.data.table(pca_obj$x[, seq_len(n_pcs_keep), drop = FALSE])
  data.table::setnames(score_values, score_cols)
  score_dt <- cbind(score_dt, score_values)
  if (n_pcs_keep < n_pcs_pad) {
    for (i in seq.int(n_pcs_keep + 1L, n_pcs_pad)) {
      score_dt[, (sprintf("PC_%d", i)) := NA_real_]
    }
  }
  score_dt
}

write_pc_summary_files <- function(pca_obj, n_pcs_keep, out_dir) {
  variance_dt <- data.table::data.table(
    pc = sprintf("PC_%d", seq_len(n_pcs_keep)),
    sdev = pca_obj$sdev[seq_len(n_pcs_keep)],
    variance_explained = (pca_obj$sdev[seq_len(n_pcs_keep)] ^ 2) / sum(pca_obj$sdev ^ 2),
    cumulative_variance = cumsum((pca_obj$sdev[seq_len(n_pcs_keep)] ^ 2) / sum(pca_obj$sdev ^ 2))
  )
  data.table::fwrite(variance_dt, file.path(out_dir, "pca_variance.csv"))

  loadings_dt <- data.table::as.data.table(pca_obj$rotation[, seq_len(n_pcs_keep), drop = FALSE], keep.rownames = "feature")
  data.table::setnames(loadings_dt, c("feature", sprintf("PC_%d", seq_len(n_pcs_keep))))
  data.table::fwrite(loadings_dt, file.path(out_dir, "pca_loadings.csv"))
}

write_representative_crops <- function(score_dt, pc_name, out_dir, n_each, n_cores = 1L) {
  required_cols <- c(
    "image_id", "filename", "label", "raw_relpath", "mask_relpath",
    "bbox_min_row", "bbox_min_col", "bbox_max_row", "bbox_max_col", pc_name
  )
  if (!all(required_cols %in% names(score_dt))) {
    missing_cols <- setdiff(required_cols, names(score_dt))
    emit_warning(sprintf(
      "Skipping crop export for %s in %s because required columns are missing: %s",
      pc_name,
      out_dir,
      paste(missing_cols, collapse = ", ")
    ))
    return(invisible(NULL))
  }

  valid_dt <- score_dt[
    is.finite(get(pc_name)) &
      !is.na(raw_relpath) & nzchar(raw_relpath) &
      !is.na(mask_relpath) & nzchar(mask_relpath)
  ]
  if (!nrow(valid_dt)) {
    emit_warning(sprintf(
      "No valid rows available for crop export on %s in %s. Rows require finite scores and non-empty raw_relpath/mask_relpath.",
      pc_name,
      out_dir
    ))
    return(invisible(NULL))
  }

  low_dt <- valid_dt[order(get(pc_name), passage_number, image_id, label)][seq_len(min(n_each, .N))]
  high_dt <- valid_dt[order(-get(pc_name), passage_number, image_id, label)][seq_len(min(n_each, .N))]

  export_side <- function(side_dt, side_name) {
    if (!nrow(side_dt)) {
      return(NULL)
    }

    side_dt[, rank_within_side := seq_len(.N)]
    by_image <- split(seq_len(nrow(side_dt)), side_dt$image_id)

    process_group <- function(idx) {
      group_dt <- side_dt[idx]
      image_path <- group_dt$raw_relpath[[1]]
      mask_path <- group_dt$mask_relpath[[1]]
      if (!file.exists(image_path) || !file.exists(mask_path)) {
        emit_warning(sprintf(
          "Missing raw or mask TIFF for image %s while exporting %s/%s. raw_relpath=%s mask_relpath=%s",
          group_dt$image_id[[1]],
          basename(dirname(out_dir)),
          side_name,
          image_path,
          mask_path
        ))
        return(NULL)
      }

      image <- read_grayscale_image(image_path)
      mask <- read_mask_tiff(mask_path)
      lapply(seq_len(nrow(group_dt)), function(j) {
        row <- group_dt[j]
        crop <- extract_object_crop(
          image = image,
          mask = mask,
          label = row$label[[1]],
          min_row = row$bbox_min_row[[1]],
          min_col = row$bbox_min_col[[1]],
          max_row = row$bbox_max_row[[1]],
          max_col = row$bbox_max_col[[1]]
        )
        list(
          crop = crop,
          meta = data.table::data.table(
            side = side_name,
            rank = as.integer(row$rank_within_side[[1]]),
            passage_number = as.integer(row$passage_number[[1]]),
            image_id = as.character(row$image_id[[1]]),
            filename = as.character(row$filename[[1]]),
            label = as.character(row$label[[1]]),
            score = as.numeric(row[[pc_name]][[1]])
          )
        )
      })
    }

    group_results <- if (.Platform$OS.type == "unix" && n_cores > 1L) {
      parallel::mclapply(by_image, process_group, mc.cores = min(n_cores, length(by_image)))
    } else {
      lapply(by_image, process_group)
    }

    flat_results <- unlist(group_results, recursive = FALSE, use.names = FALSE)
    flat_results <- Filter(Negate(is.null), flat_results)
    if (!length(flat_results)) {
      emit_warning(sprintf("No crops could be extracted for %s in %s/%s.", pc_name, out_dir, side_name))
      return(NULL)
    }

    crop_mats <- lapply(flat_results, `[[`, "crop")
    meta_dt <- data.table::rbindlist(lapply(flat_results, `[[`, "meta"), use.names = TRUE, fill = TRUE)
    panel <- assemble_crop_panel(crop_mats)
    panel_path <- file.path(out_dir, sprintf("%s_panel.tiff", side_name))
    write_crop_tiff(panel, panel_path)
    data.table::fwrite(meta_dt, file.path(out_dir, sprintf("%s_panel_metadata.csv", side_name)))
    meta_dt[, panel_path := panel_path]
    meta_dt
  }

  low_meta <- export_side(low_dt, "low")
  high_meta <- export_side(high_dt, "high")
  combined_meta <- data.table::rbindlist(list(low_meta, high_meta), use.names = TRUE, fill = TRUE)
  if (nrow(combined_meta)) {
    data.table::fwrite(combined_meta, file.path(out_dir, "panel_metadata.csv"))
  }
}

main <- function() {
  required_packages <- c("data.table", "ggplot2", "tiff", "parallel")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      sprintf("Install required packages before running: %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }

  args <- parse_args(commandArgs(trailingOnly = TRUE))
  project_root <- resolve_project_root()
  if (!is.finite(args$n_cores) || args$n_cores < 1L) {
    args$n_cores <- 1L
  }

  lineage_rds <- project_path(project_root, args$lineage_rds)
  segmentation_manifest <- project_path(project_root, args$segmentation_manifest)
  embedding_dir <- project_path(project_root, args$embedding_dir)
  out_dir <- file.path(project_path(project_root, args$out_dir), paste0("lineage_", args$lineage_id))
  ensure_dir(out_dir)
  reset_warning_state(file.path(out_dir, "warning_log.txt"))

  if (!file.exists(lineage_rds)) {
    stop(sprintf("Lineage RDS not found: %s", lineage_rds), call. = FALSE)
  }
  if (!file.exists(segmentation_manifest)) {
    stop(sprintf("Segmentation manifest not found: %s", segmentation_manifest), call. = FALSE)
  }

  lineages <- readRDS(lineage_rds)
  lineage_obj <- extract_lineage(lineages, args$lineage_id)
  image_meta <- build_lineage_image_meta(lineage_obj)

  seg_dt <- data.table::fread(segmentation_manifest)
  if ("id" %in% names(seg_dt)) {
    seg_dt[, id := as.character(id)]
  }
  if ("filename" %in% names(seg_dt)) {
    seg_dt[, filename := as.character(filename)]
  }
  if ("source_filepath" %in% names(seg_dt)) {
    seg_dt[, source_filepath := as.character(source_filepath)]
    seg_dt[, source_basename := basename(source_filepath)]
  }
  if ("source_filepath" %in% names(seg_dt)) {
    seg_dt[, source_filepath := project_path(project_root, normalize_relpath(source_filepath))]
  }
  for (path_col in intersect(c("raw_relpath", "mask_relpath", "embedding_relpath"), names(seg_dt))) {
    seg_dt[, (path_col) := project_path(project_root, normalize_relpath(get(path_col)))]
  }

  resolved_meta_list <- lapply(seq_len(nrow(image_meta)), function(i) {
    row <- image_meta[i]
    seg_row <- match_segmentation_row(seg_dt, row)

    embedding_path <- NA_character_
    raw_path <- NA_character_
    mask_path <- NA_character_

    if (!is.null(seg_row) && nrow(seg_row)) {
      if ("embedding_relpath" %in% names(seg_row) && file.exists(seg_row$embedding_relpath[[1]])) {
        embedding_path <- seg_row$embedding_relpath[[1]]
      }
      if ("raw_relpath" %in% names(seg_row) && file.exists(seg_row$raw_relpath[[1]])) {
        raw_path <- seg_row$raw_relpath[[1]]
      } else if ("source_filepath" %in% names(seg_row) && file.exists(seg_row$source_filepath[[1]])) {
        raw_path <- seg_row$source_filepath[[1]]
      }
      if ("mask_relpath" %in% names(seg_row) && file.exists(seg_row$mask_relpath[[1]])) {
        mask_path <- seg_row$mask_relpath[[1]]
      }
    }

    if (!isTRUE(nzchar(embedding_path)) || is.na(embedding_path)) {
      candidates <- candidate_embedding_paths(
        project_root = project_root,
        embedding_dir = embedding_dir,
        image_id = row$image_id[[1]],
        filename = if ("filename" %in% names(row)) row$filename[[1]] else NA_character_
      )
      existing <- candidates[file.exists(candidates)]
      if (length(existing)) {
        embedding_path <- existing[[1]]
      }
    }

    data.table::data.table(
      image_id = as.character(row$image_id[[1]]),
      filename = if ("filename" %in% names(row)) as.character(row$filename[[1]]) else NA_character_,
      filepath = if ("filepath" %in% names(row)) as.character(row$filepath[[1]]) else NA_character_,
      passage_id = as.character(row$passage_id[[1]]),
      passage_number = as.integer(row$passage_number[[1]]),
      embedding_path = embedding_path,
      raw_relpath = raw_path,
      mask_relpath = mask_path
    )
  })
  resolved_meta <- data.table::rbindlist(resolved_meta_list, use.names = TRUE, fill = TRUE)
  resolved_meta_all <- data.table::copy(resolved_meta)
  duplicate_meta <- resolved_meta_all[, .N, by = .(image_id)][N > 1L]
  if (nrow(duplicate_meta)) {
    emit_warning(sprintf(
      "Collapsing duplicated lineage image metadata for %d image_id values: %s",
      nrow(duplicate_meta),
      paste(duplicate_meta$image_id, collapse = ", ")
    ))
  }
  resolved_meta <- resolved_meta_all[order(image_id, passage_number)][, .SD[1], by = image_id]

  missing_embedding <- resolved_meta[is.na(embedding_path) | !nzchar(embedding_path)]
  if (nrow(missing_embedding)) {
    emit_warning(sprintf(
      "Missing embedding CSV for %d lineage images: %s",
      nrow(missing_embedding),
      paste(missing_embedding$image_id, collapse = ", ")
    ))
  }
  missing_raw_mask <- resolved_meta[is.na(raw_relpath) | !nzchar(raw_relpath) | is.na(mask_relpath) | !nzchar(mask_relpath)]
  if (nrow(missing_raw_mask)) {
    emit_warning(sprintf(
      "Missing raw/mask paths for %d lineage images after segmentation-manifest matching: %s",
      nrow(missing_raw_mask),
      paste(missing_raw_mask$image_id, collapse = ", ")
    ))
  }
  cat(sprintf(
    "Resolved image paths: embeddings=%d/%d raw_mask=%d/%d\n",
    sum(!is.na(resolved_meta$embedding_path) & nzchar(resolved_meta$embedding_path)),
    nrow(resolved_meta),
    sum(!is.na(resolved_meta$raw_relpath) & nzchar(resolved_meta$raw_relpath) & !is.na(resolved_meta$mask_relpath) & nzchar(resolved_meta$mask_relpath)),
    nrow(resolved_meta)
  ))

  embedding_tables <- lapply(seq_len(nrow(resolved_meta)), function(i) {
    row <- resolved_meta[i]
    if (is.na(row$embedding_path[[1]]) || !nzchar(row$embedding_path[[1]]) || !file.exists(row$embedding_path[[1]])) {
      return(NULL)
    }
    dt <- data.table::fread(row$embedding_path[[1]])
    if (!nrow(dt)) {
      dt <- data.table::data.table()
    }
    if (nrow(dt)) {
      dt[, image_id := row$image_id[[1]]]
      dt[, filename := row$filename[[1]]]
      dt[, passage_id := row$passage_id[[1]]]
      dt[, passage_number := row$passage_number[[1]]]
      dt[, raw_relpath := row$raw_relpath[[1]]]
      dt[, mask_relpath := row$mask_relpath[[1]]]
    }
    dt
  })
  object_dt <- data.table::rbindlist(embedding_tables, use.names = TRUE, fill = TRUE)

  if (!nrow(object_dt)) {
    stop(sprintf("No embedding rows could be loaded for lineage %s.", args$lineage_id), call. = FALSE)
  }

  image_counts <- merge(
    resolved_meta[, .(image_id, filename, passage_id, passage_number, embedding_path, raw_relpath, mask_relpath)],
    object_dt[, .(cell_count = .N), by = .(image_id)],
    by = "image_id",
    all.x = TRUE
  )
  image_counts[!file.exists(embedding_path) | is.na(embedding_path) | !nzchar(embedding_path), cell_count := NA_integer_]
  image_counts[is.na(cell_count) & !is.na(embedding_path) & nzchar(embedding_path), cell_count := 0L]

  available_counts <- image_counts[is.finite(cell_count), cell_count]
  if (!length(available_counts)) {
    stop(sprintf("No valid image cell counts available for lineage %s.", args$lineage_id), call. = FALSE)
  }
  cutpoints <- as.numeric(stats::quantile(available_counts, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE, type = 7))
  image_counts[, quantile := vapply(cell_count, assign_count_quantile, character(1), cutpoints = cutpoints)]
  image_counts[, quantile := factor(quantile, levels = c("Q1_low", "Q2_mid", "Q3_high"))]
  data.table::fwrite(image_counts, file.path(out_dir, "image_cell_counts.csv"))
  build_image_count_plot(
    image_counts = image_counts,
    cutpoints = cutpoints,
    out_path = file.path(out_dir, "image_cell_count_distribution.png"),
    lineage_id = args$lineage_id
  )

  image_quantiles <- unique(
    image_counts[, .(image_id, quantile, cell_count)],
    by = "image_id"
  )

  object_dt <- merge(
    object_dt,
    image_quantiles,
    by = "image_id",
    all.x = TRUE
  )

  feature_cols <- c(intersect("log1p_area", names(object_dt)), grep("^embedding_", names(object_dt), value = TRUE))
  if (!length(feature_cols)) {
    stop("No PCA feature columns found. Expected `log1p_area` and/or `embedding_*` columns.", call. = FALSE)
  }

  combined_score_tables <- list()
  quantile_levels <- c("Q1_low", "Q2_mid", "Q3_high")

  for (quantile_name in quantile_levels) {
    quantile_out_dir <- file.path(out_dir, quantile_name)
    ensure_dir(quantile_out_dir)

    quantile_dt <- object_dt[quantile == quantile_name]
    if (!nrow(quantile_dt)) {
      emit_warning(sprintf("No cell rows available for %s in lineage %s.", quantile_name, args$lineage_id))
      next
    }

    feature_mat <- as.matrix(quantile_dt[, ..feature_cols])
    storage.mode(feature_mat) <- "double"
    keep_rows <- stats::complete.cases(feature_mat)
    quantile_complete <- quantile_dt[keep_rows]
    feature_mat <- feature_mat[keep_rows, , drop = FALSE]

    if (nrow(feature_mat) < 3L) {
      emit_warning(sprintf("Too few complete rows for PCA in %s of lineage %s.", quantile_name, args$lineage_id))
      next
    }

    pca_obj <- stats::prcomp(feature_mat, center = TRUE, scale. = TRUE)
    n_pcs_keep <- min(args$n_pcs, ncol(pca_obj$x))
    score_dt <- build_score_table(
      object_dt = quantile_complete,
      pca_obj = pca_obj,
      n_pcs_keep = n_pcs_keep,
      n_pcs_pad = args$n_pcs
    )
    score_dt[, quantile := quantile_name]

    write_pc_summary_files(pca_obj, n_pcs_keep, quantile_out_dir)
    saveRDS(pca_obj, file.path(quantile_out_dir, "pca_model.rds"))

    cor_dt <- data.table::data.table(
      pc = sprintf("PC_%d", seq_len(n_pcs_keep)),
      correlation = vapply(seq_len(n_pcs_keep), function(i) {
        safe_cor(score_dt[[sprintf("PC_%d", i)]], score_dt$passage_number)
      }, numeric(1))
    )
    cor_dt[, abs_correlation := abs(correlation)]
    cor_dt <- cor_dt[order(-abs_correlation, pc)]
    data.table::fwrite(cor_dt, file.path(quantile_out_dir, "pc_passage_correlations.csv"))

    top_corr_pcs <- head(cor_dt$pc, min(args$top_corr_pcs, nrow(cor_dt)))
    top_export_pcs <- head(cor_dt$pc, min(args$top_export_pcs, nrow(cor_dt)))

    for (pc_name in top_corr_pcs) {
        write_representative_crops(
          score_dt = score_dt,
          pc_name = pc_name,
          out_dir = file.path(quantile_out_dir, "top3_abs_correlation", pc_name),
          n_each = args$crops_per_side,
          n_cores = args$n_cores
        )
      }

    for (pc_name in top_export_pcs) {
        write_representative_crops(
          score_dt = score_dt,
          pc_name = pc_name,
          out_dir = file.path(quantile_out_dir, "top5_abs_correlation", pc_name),
          n_each = args$crops_per_side,
          n_cores = args$n_cores
        )
      }

    keep_cols <- c("image_id", "quantile", "label", sprintf("PC_%d", seq_len(args$n_pcs)))
    extra_cols <- intersect(c("filename", "passage_id", "passage_number"), names(score_dt))
    select_cols <- c(keep_cols, extra_cols)
    score_export <- score_dt[, ..select_cols]
    data.table::fwrite(score_export, file.path(quantile_out_dir, "pc_scores.csv"))
    combined_score_tables[[quantile_name]] <- score_export
  }

  if (!length(combined_score_tables)) {
    stop(sprintf("No quantile produced PCA outputs for lineage %s.", args$lineage_id), call. = FALSE)
  }

  combined_scores <- data.table::rbindlist(combined_score_tables, use.names = TRUE, fill = TRUE)
  data.table::fwrite(combined_scores, file.path(out_dir, "lineage_pc_scores.csv"))
  saveRDS(combined_scores, file.path(out_dir, "lineage_pc_scores.rds"))

  cat(sprintf("Lineage: %s\n", args$lineage_id))
  cat(sprintf("Images in lineage: %d\n", nrow(image_counts)))
  cat(sprintf("Cells loaded: %d\n", nrow(object_dt)))
  cat(sprintf("Crop export cores: %d\n", args$n_cores))
  cat(sprintf("Warnings emitted: %d\n", length(warning_state$messages)))
  if (length(warning_state$messages)) {
    cat(sprintf("Warning log: %s\n", normalize_relpath(warning_state$log_path)))
  }
  cat(sprintf("Output written to %s\n", normalize_relpath(out_dir)))
}

main()
