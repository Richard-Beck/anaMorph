parse_args <- function(args) {
  parsed <- list(
    lineage_rds = "core_data/lineages.Rds",
    lineage_scores_dir = "data/lineage_embeddings",
    segmentation_manifest = "data/all_images/manifests/segmentation_images.csv",
    embedding_dir = "data/all_images/embeddings",
    out_dir = "data/lineage_embedding_clusters",
    k_max = "20",
    nstart = "20",
    representative_cells = "12",
    n_cores = "1",
    max_metric_cells = "5000"
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

  parsed$k_max <- as.integer(parsed$k_max)
  parsed$nstart <- as.integer(parsed$nstart)
  parsed$representative_cells <- as.integer(parsed$representative_cells)
  parsed$n_cores <- as.integer(parsed$n_cores)
  parsed$max_metric_cells <- as.integer(parsed$max_metric_cells)
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

match_segmentation_row <- function(seg_dt, image_row) {
  seg_row <- NULL

  if ("filename" %in% names(seg_dt) && "filename" %in% names(image_row)) {
    seg_row <- seg_dt[filename == image_row$filename[[1]]]
    if (nrow(seg_row)) return(seg_row[1])
  }

  if ("source_filepath" %in% names(seg_dt) && "filepath" %in% names(image_row)) {
    seg_row <- seg_dt[source_filepath == image_row$filepath[[1]]]
    if (nrow(seg_row)) return(seg_row[1])
  }

  if ("source_basename" %in% names(seg_dt) && "filepath" %in% names(image_row)) {
    seg_row <- seg_dt[source_basename == basename(image_row$filepath[[1]])]
    if (nrow(seg_row)) return(seg_row[1])
  }

  if ("id" %in% names(seg_dt) && "sample_id" %in% names(image_row)) {
    seg_row <- seg_dt[id == image_row$sample_id[[1]]]
    if (nrow(seg_row)) return(seg_row[1])
  }

  NULL
}

candidate_embedding_paths <- function(embedding_dir, image_id, filename = NA_character_) {
  stems <- unique(stats::na.omit(c(
    as.character(image_id),
    tools::file_path_sans_ext(basename(as.character(filename)))
  )))
  file.path(embedding_dir, paste0(stems, "_embeddings.csv"))
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
  if (!is.finite(hi) || hi <= lo) hi <- max(x, na.rm = TRUE)
  if (!is.finite(hi) || hi <= lo) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x)))
  }
  out <- (x - lo) / (hi - lo)
  out[out < 0] <- 0
  out[out > 1] <- 1
  out
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
  if (!length(crop_mats)) return(NULL)
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

evaluate_k_grid <- function(dt, feature_cols, k_max, nstart, seed = 1L, max_metric_cells = 5000L) {
  feature_mat <- as.matrix(dt[, ..feature_cols])
  storage.mode(feature_mat) <- "double"
  keep_rows <- stats::complete.cases(feature_mat)
  filtered_dt <- dt[keep_rows]
  filtered_mat <- feature_mat[keep_rows, , drop = FALSE]
  n_obs <- nrow(filtered_mat)

  if (n_obs < 3L) {
    stop("Need at least three complete rows for clustering.", call. = FALSE)
  }

  candidate_ks <- seq.int(2L, min(as.integer(k_max), n_obs - 1L))
  if (!length(candidate_ks)) {
    stop("No valid cluster counts available.", call. = FALSE)
  }

  metric_idx <- seq_len(n_obs)
  if (n_obs > max_metric_cells) {
    set.seed(seed)
    metric_idx <- sort(sample.int(n_obs, size = max_metric_cells, replace = FALSE))
  }
  metric_mat <- filtered_mat[metric_idx, , drop = FALSE]
  metric_n <- nrow(metric_mat)
  dist_mat <- stats::dist(metric_mat)
  eval_list <- lapply(candidate_ks, function(k) {
    set.seed(seed)
    km <- stats::kmeans(metric_mat, centers = k, nstart = nstart)
    sil <- cluster::silhouette(km$cluster, dist_mat)
    mean_sil <- mean(sil[, "sil_width"])
    ch <- (km$betweenss / (k - 1)) / (km$tot.withinss / (metric_n - k))
    list(
      k = k,
      model = km,
      mean_silhouette = mean_sil,
      calinski_harabasz = ch,
      tot_withinss = km$tot.withinss
    )
  })

  metrics_dt <- data.table::rbindlist(lapply(eval_list, function(x) {
    data.table::data.table(
      k = x$k,
      mean_silhouette = x$mean_silhouette,
      calinski_harabasz = x$calinski_harabasz,
      tot_withinss = x$tot_withinss
    )
  }))

  best_idx <- which.max(metrics_dt$mean_silhouette)
  best_result <- eval_list[[best_idx]]
  set.seed(seed)
  final_model <- stats::kmeans(filtered_mat, centers = best_result$k, nstart = nstart)
  centers <- final_model$centers[final_model$cluster, , drop = FALSE]
  cluster_distance <- sqrt(rowSums((filtered_mat - centers) ^ 2))

  clustered_dt <- data.table::copy(filtered_dt)
  clustered_dt[, cluster := factor(final_model$cluster)]
  clustered_dt[, cluster_distance := cluster_distance]

  list(
    metrics = metrics_dt,
    clustered_tbl = clustered_dt,
    best_k = best_result$k,
    best_model = final_model,
    feature_cols = feature_cols,
    metric_cells = metric_n,
    clustered_cells = n_obs
  )
}

plot_metric_grid <- function(metrics_dt, title) {
  plot_dt <- data.table::melt(
    data.table::copy(metrics_dt),
    id.vars = "k",
    measure.vars = c("mean_silhouette", "calinski_harabasz", "tot_withinss"),
    variable.name = "metric",
    value.name = "value"
  )

  ggplot2::ggplot(plot_dt, ggplot2::aes(x = k, y = value)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = "Number of clusters", y = "Metric value")
}

build_cluster_frequency <- function(clustered_dt) {
  freq_dt <- clustered_dt[, .N, by = .(passage_number, cluster)]
  total_dt <- clustered_dt[, .(total = .N), by = .(passage_number)]
  freq_dt <- merge(freq_dt, total_dt, by = "passage_number", all.x = TRUE)
  freq_dt[, relative_frequency := N / total]
  freq_dt
}

plot_cluster_frequency <- function(freq_dt, title) {
  ggplot2::ggplot(freq_dt, ggplot2::aes(x = passage_number, y = relative_frequency, color = cluster)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = "Passage number", y = "Relative cluster frequency")
}

cluster_frequency_correlations <- function(freq_dt) {
  freq_dt[, .(
    correlation = if (.N >= 3L) stats::cor(passage_number, relative_frequency) else NA_real_
  ), by = cluster][order(cluster)]
}

resolve_image_metadata <- function(score_dt, lineage_obj, seg_dt, embedding_dir) {
  lineage_meta <- data.table::rbindlist(lineage_obj$passages, use.names = TRUE, fill = TRUE)
  lineage_meta <- data.table::as.data.table(lineage_meta)
  if ("id" %in% names(lineage_meta)) lineage_meta[, id := as.character(id)]
  if ("filename" %in% names(lineage_meta)) lineage_meta[, filename := as.character(filename)]
  if ("filepath" %in% names(lineage_meta)) lineage_meta[, filepath := as.character(filepath)]
  lineage_meta[, image_id := if ("filename" %in% names(lineage_meta)) filename else id]

  score_images <- unique(score_dt[, .(image_id, filename)])
  score_images <- merge(score_images, lineage_meta[, .(image_id, sample_id = id, filepath)], by = "image_id", all.x = TRUE)

  resolved_list <- lapply(seq_len(nrow(score_images)), function(i) {
    row <- score_images[i]
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
      candidates <- candidate_embedding_paths(embedding_dir, row$image_id[[1]], row$filename[[1]])
      existing <- candidates[file.exists(candidates)]
      if (length(existing)) embedding_path <- existing[[1]]
    }

    data.table::data.table(
      image_id = as.character(row$image_id[[1]]),
      filename = as.character(row$filename[[1]]),
      sample_id = as.character(row$sample_id[[1]]),
      filepath = as.character(row$filepath[[1]]),
      embedding_path = embedding_path,
      raw_relpath = raw_path,
      mask_relpath = mask_path
    )
  })

  data.table::rbindlist(resolved_list, use.names = TRUE, fill = TRUE)
}

load_crop_metadata <- function(score_dt, resolved_meta) {
  crop_sources <- resolved_meta[!is.na(embedding_path) & nzchar(embedding_path)]
  crop_tables <- lapply(seq_len(nrow(crop_sources)), function(i) {
    row <- crop_sources[i]
    dt <- data.table::fread(row$embedding_path[[1]], select = c("label", "bbox_min_row", "bbox_min_col", "bbox_max_row", "bbox_max_col"))
    dt[, image_id := row$image_id[[1]]]
    dt[, raw_relpath := row$raw_relpath[[1]]]
    dt[, mask_relpath := row$mask_relpath[[1]]]
    dt
  })
  crop_dt <- data.table::rbindlist(crop_tables, use.names = TRUE, fill = TRUE)
  crop_dt[, label := as.character(label)]
  score_key_dt <- unique(score_dt[, .(image_id, label = as.character(label))])
  merge(score_key_dt, crop_dt, by = c("image_id", "label"), all.x = TRUE)
}

write_cluster_panels <- function(clustered_dt, crop_meta_dt, out_dir, n_each, n_cores = 1L) {
  if (!nrow(clustered_dt)) return(invisible(NULL))
  merged_dt <- merge(
    clustered_dt,
    crop_meta_dt,
    by = c("image_id", "label"),
    all.x = TRUE
  )

  missing_crop <- merged_dt[is.na(raw_relpath) | !nzchar(raw_relpath) | is.na(mask_relpath) | !nzchar(mask_relpath)]
  if (nrow(missing_crop)) {
    emit_warning(sprintf(
      "Missing crop metadata for %d clustered rows in %s.",
      nrow(missing_crop),
      out_dir
    ))
  }

  cluster_levels <- sort(unique(as.character(merged_dt$cluster)))
  summary_rows <- list()

  for (cluster_id in cluster_levels) {
    cluster_dt <- merged_dt[cluster == cluster_id & !is.na(raw_relpath) & nzchar(raw_relpath) & !is.na(mask_relpath) & nzchar(mask_relpath)]
    if (!nrow(cluster_dt)) {
      emit_warning(sprintf("No valid crop rows for cluster %s in %s.", cluster_id, out_dir))
      next
    }

    cluster_dt <- cluster_dt[order(cluster_distance, passage_number, image_id, label)][seq_len(min(n_each, .N))]
    by_image <- split(seq_len(nrow(cluster_dt)), cluster_dt$image_id)

    process_group <- function(idx) {
      group_dt <- cluster_dt[idx]
      image_path <- group_dt$raw_relpath[[1]]
      mask_path <- group_dt$mask_relpath[[1]]
      if (!file.exists(image_path) || !file.exists(mask_path)) {
        emit_warning(sprintf(
          "Missing raw or mask TIFF for image %s while exporting cluster %s. raw_relpath=%s mask_relpath=%s",
          group_dt$image_id[[1]],
          cluster_id,
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
            cluster = as.character(cluster_id),
            rank = j,
            passage_number = as.integer(row$passage_number[[1]]),
            image_id = as.character(row$image_id[[1]]),
            filename = as.character(row$filename[[1]]),
            label = as.character(row$label[[1]]),
            cluster_distance = as.numeric(row$cluster_distance[[1]])
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
      emit_warning(sprintf("No crops extracted for cluster %s in %s.", cluster_id, out_dir))
      next
    }

    crop_mats <- lapply(flat_results, `[[`, "crop")
    panel <- assemble_crop_panel(crop_mats)
    cluster_dir <- file.path(out_dir, sprintf("cluster_%s", cluster_id))
    ensure_dir(cluster_dir)
    panel_path <- file.path(cluster_dir, "representative_cells.tiff")
    write_crop_tiff(panel, panel_path)

    meta_dt <- data.table::rbindlist(lapply(flat_results, `[[`, "meta"), use.names = TRUE, fill = TRUE)
    meta_dt[, panel_path := panel_path]
    data.table::fwrite(meta_dt, file.path(cluster_dir, "representative_cells.csv"))
    summary_rows[[cluster_id]] <- meta_dt
  }

  summary_dt <- data.table::rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
  if (nrow(summary_dt)) {
    data.table::fwrite(summary_dt, file.path(out_dir, "cluster_representative_cells.csv"))
  }
}

run_scope_clustering <- function(scope_name, score_dt, feature_cols, crop_meta_dt, out_dir, k_max, nstart, representative_cells, n_cores, max_metric_cells) {
  ensure_dir(out_dir)
  result <- evaluate_k_grid(
    score_dt,
    feature_cols,
    k_max = k_max,
    nstart = nstart,
    max_metric_cells = max_metric_cells
  )

  data.table::fwrite(result$metrics, file.path(out_dir, "clustering_metrics.csv"))
  metric_plot <- plot_metric_grid(result$metrics, sprintf("%s clustering metrics", scope_name))
  ggplot2::ggsave(file.path(out_dir, "clustering_metrics.png"), metric_plot, width = 8, height = 8, dpi = 150)

  clustered_dt <- data.table::copy(result$clustered_tbl)
  data.table::fwrite(clustered_dt, file.path(out_dir, "best_clustering.csv"))

  freq_dt <- build_cluster_frequency(clustered_dt)
  cor_dt <- cluster_frequency_correlations(freq_dt)
  data.table::fwrite(freq_dt, file.path(out_dir, "cluster_frequency_by_passage.csv"))
  data.table::fwrite(cor_dt, file.path(out_dir, "cluster_frequency_correlations.csv"))

  freq_plot <- plot_cluster_frequency(freq_dt, sprintf("%s cluster frequencies over time", scope_name))
  ggplot2::ggsave(file.path(out_dir, "cluster_frequency_over_time.png"), freq_plot, width = 9, height = 5.5, dpi = 150)

  write_cluster_panels(
    clustered_dt = clustered_dt,
    crop_meta_dt = crop_meta_dt,
    out_dir = file.path(out_dir, "representative_clusters"),
    n_each = representative_cells,
    n_cores = n_cores
  )

  list(
    best_k = result$best_k,
    n_cells = nrow(clustered_dt),
    metrics = result$metrics,
    metric_cells = result$metric_cells
  )
}

main <- function() {
  required_packages <- c("data.table", "ggplot2", "cluster", "parallel", "tiff")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      sprintf("Install required packages before running: %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }

  args <- parse_args(commandArgs(trailingOnly = TRUE))
  project_root <- resolve_project_root()
  if (!is.finite(args$n_cores) || args$n_cores < 1L) args$n_cores <- 1L
  if (!is.finite(args$max_metric_cells) || args$max_metric_cells < 100L) args$max_metric_cells <- 5000L

  lineage_rds <- project_path(project_root, args$lineage_rds)
  lineage_scores_dir <- file.path(project_path(project_root, args$lineage_scores_dir), paste0("lineage_", args$lineage_id))
  segmentation_manifest <- project_path(project_root, args$segmentation_manifest)
  embedding_dir <- project_path(project_root, args$embedding_dir)
  out_dir <- file.path(project_path(project_root, args$out_dir), paste0("lineage_", args$lineage_id))
  ensure_dir(out_dir)
  reset_warning_state(file.path(out_dir, "warning_log.txt"))

  lineages <- readRDS(lineage_rds)
  lineage_obj <- lineages[[args$lineage_id]]
  if (is.null(lineage_obj)) {
    hit <- vapply(lineages, function(x) identical(as.character(x$lineage_id), as.character(args$lineage_id)), logical(1))
    if (!any(hit)) stop(sprintf("Lineage not found: %s", args$lineage_id), call. = FALSE)
    lineage_obj <- lineages[[which(hit)[[1]]]]
  }

  score_path <- file.path(lineage_scores_dir, "lineage_pc_scores.csv")
  if (!file.exists(score_path)) {
    stop(sprintf("Lineage PC score file not found: %s", score_path), call. = FALSE)
  }

  score_dt <- data.table::fread(score_path)
  if (!nrow(score_dt)) {
    stop("No PC scores found for clustering.", call. = FALSE)
  }
  score_dt[, image_id := as.character(image_id)]
  score_dt[, label := as.character(label)]
  if ("filename" %in% names(score_dt)) score_dt[, filename := as.character(filename)]
  if ("quantile" %in% names(score_dt)) score_dt[, quantile := as.character(quantile)]

  feature_cols <- grep("^PC_[0-9]+$", names(score_dt), value = TRUE)
  if (!length(feature_cols)) {
    stop("No PC feature columns found. Expected columns like `PC_1`.", call. = FALSE)
  }

  seg_dt <- data.table::fread(segmentation_manifest)
  if ("id" %in% names(seg_dt)) seg_dt[, id := as.character(id)]
  if ("filename" %in% names(seg_dt)) seg_dt[, filename := as.character(filename)]
  if ("source_filepath" %in% names(seg_dt)) {
    seg_dt[, source_filepath := as.character(source_filepath)]
    seg_dt[, source_basename := basename(source_filepath)]
    seg_dt[, source_filepath := project_path(project_root, normalize_relpath(source_filepath))]
  }
  for (path_col in intersect(c("raw_relpath", "mask_relpath", "embedding_relpath"), names(seg_dt))) {
    seg_dt[, (path_col) := project_path(project_root, normalize_relpath(get(path_col)))]
  }

  resolved_meta <- resolve_image_metadata(score_dt, lineage_obj, seg_dt, embedding_dir)
  cat(sprintf(
    "Resolved image paths for clustering: embeddings=%d/%d raw_mask=%d/%d\n",
    sum(!is.na(resolved_meta$embedding_path) & nzchar(resolved_meta$embedding_path)),
    nrow(resolved_meta),
    sum(!is.na(resolved_meta$raw_relpath) & nzchar(resolved_meta$raw_relpath) & !is.na(resolved_meta$mask_relpath) & nzchar(resolved_meta$mask_relpath)),
    nrow(resolved_meta)
  ))

  crop_meta_dt <- load_crop_metadata(score_dt, resolved_meta)
  cluster_summaries <- list()

  pooled_out_dir <- file.path(out_dir, "pooled_all_quantiles")
  cluster_summaries[["pooled_all_quantiles"]] <- run_scope_clustering(
    scope_name = sprintf("%s pooled across quantiles", args$lineage_id),
    score_dt = score_dt,
    feature_cols = feature_cols,
    crop_meta_dt = crop_meta_dt,
    out_dir = pooled_out_dir,
    k_max = args$k_max,
    nstart = args$nstart,
    representative_cells = args$representative_cells,
    n_cores = args$n_cores,
    max_metric_cells = args$max_metric_cells
  )

  for (quantile_name in sort(unique(stats::na.omit(score_dt$quantile)))) {
    quantile_dt <- score_dt[quantile == quantile_name]
    if (nrow(quantile_dt) < 3L) {
      emit_warning(sprintf("Skipping %s because it has fewer than three rows.", quantile_name))
      next
    }
    cluster_summaries[[quantile_name]] <- run_scope_clustering(
      scope_name = sprintf("%s %s", args$lineage_id, quantile_name),
      score_dt = quantile_dt,
      feature_cols = feature_cols,
      crop_meta_dt = crop_meta_dt,
      out_dir = file.path(out_dir, quantile_name),
      k_max = args$k_max,
      nstart = args$nstart,
      representative_cells = args$representative_cells,
      n_cores = args$n_cores,
      max_metric_cells = args$max_metric_cells
    )
  }

  summary_dt <- data.table::rbindlist(lapply(names(cluster_summaries), function(scope_name) {
    x <- cluster_summaries[[scope_name]]
    data.table::data.table(
      scope = scope_name,
      best_k = x$best_k,
      n_cells = x$n_cells,
      metric_cells = x$metric_cells
    )
  }), use.names = TRUE, fill = TRUE)
  data.table::fwrite(summary_dt, file.path(out_dir, "clustering_summary.csv"))

  cat(sprintf("Lineage: %s\n", args$lineage_id))
  cat(sprintf("Cells available for clustering: %d\n", nrow(score_dt)))
  cat(sprintf("PC features used: %d\n", length(feature_cols)))
  cat(sprintf("Metric evaluation cells (max): %d\n", args$max_metric_cells))
  cat(sprintf("Crop export cores: %d\n", args$n_cores))
  cat(sprintf("Warnings emitted: %d\n", length(warning_state$messages)))
  if (length(warning_state$messages)) {
    cat(sprintf("Warning log: %s\n", normalize_relpath(warning_state$log_path)))
  }
  cat(sprintf("Output written to %s\n", normalize_relpath(out_dir)))
}

main()
