resolve_project_root <- function() {
  candidates <- c(
    normalizePath(".", winslash = "/", mustWork = FALSE),
    normalizePath("..", winslash = "/", mustWork = FALSE),
    normalizePath("../..", winslash = "/", mustWork = FALSE)
  )

  for (candidate in unique(candidates)) {
    if (file.exists(file.path(candidate, "README.md")) &&
        dir.exists(file.path(candidate, "manifests"))) {
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

embedding_column_names <- function(x) {
  grep("^embedding_[0-9]{3}$", names(x), value = TRUE)
}

load_embedding_vectors <- function(
  image_names,
  segmentation_manifest = "manifests/dev_subset_segmentation_images.csv",
  metadata_manifest = "manifests/dev_subset_images.csv"
) {
  required_packages <- c("data.table")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      sprintf("Install required packages before running: %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }

  if (length(image_names) == 0) {
    stop("`image_names` must contain at least one filename.", call. = FALSE)
  }

  image_names <- unique(as.character(image_names))
  project_root <- resolve_project_root()
  segmentation_manifest <- project_path(project_root, segmentation_manifest)
  metadata_manifest <- project_path(project_root, metadata_manifest)

  metadata_dt <- data.table::fread(metadata_manifest, select = c("filename", "id"))
  requested_dt <- metadata_dt[filename %in% image_names]

  unknown_images <- setdiff(image_names, requested_dt$filename)
  if (length(unknown_images) > 0) {
    stop(
      sprintf(
        "No image metadata entries found for: %s",
        paste(unknown_images, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  requested_ids <- unique(requested_dt$id)

  seg_dt <- data.table::fread(segmentation_manifest)
  seg_dt[, embedding_relpath := project_path(project_root, normalize_relpath(embedding_relpath))]

  resolved_dt <- unique(seg_dt[id %in% requested_ids, .(
    id,
    filename,
    embedding_relpath
  )], by = "filename")

  missing_ids <- setdiff(requested_ids, unique(resolved_dt$id))
  if (length(missing_ids) > 0) {
    stop(
      sprintf(
        "No embedding manifest entries found for id values: %s",
        paste(missing_ids, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  tables <- lapply(seq_len(nrow(resolved_dt)), function(i) {
    seg_row <- resolved_dt[i]
    dt <- data.table::fread(seg_row$embedding_relpath)
    emb_cols <- embedding_column_names(dt)
    if (length(emb_cols) == 0) {
      stop(sprintf("No embedding columns found in %s", seg_row$embedding_relpath), call. = FALSE)
    }

    out <- dt[, c("label", emb_cols), with = FALSE]
    out[, image_name := seg_row$filename]
    out[, id := seg_row$id]
    data.table::setcolorder(out, c("image_name", "id", "label", emb_cols))
    out
  })

  data.table::rbindlist(tables, use.names = TRUE, fill = FALSE)
}

run_embedding_pca <- function(embedding_tbl, embedding_cols = NULL, n_pcs = 20, center = TRUE, scale. = TRUE) {
  required_packages <- c("data.table")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      sprintf("Install required packages before running: %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }

  dt <- data.table::as.data.table(embedding_tbl)
  if (is.null(embedding_cols)) {
    embedding_cols <- embedding_column_names(dt)
  }

  if (length(embedding_cols) == 0) {
    stop("No embedding columns found. Expected columns matching `^embedding_[0-9]{3}$`.", call. = FALSE)
  }

  feature_dt <- dt[, ..embedding_cols]
  feature_mat <- as.matrix(feature_dt)
  storage.mode(feature_mat) <- "double"

  keep_rows <- stats::complete.cases(feature_mat)
  if (!any(keep_rows)) {
    stop("No complete embedding rows available for PCA.", call. = FALSE)
  }

  filtered_mat <- feature_mat[keep_rows, , drop = FALSE]
  filtered_dt <- dt[keep_rows]
  pca <- stats::prcomp(filtered_mat, center = center, scale. = scale.)
  n_pcs <- min(as.integer(n_pcs), ncol(pca$x))
  if (n_pcs < 1) {
    stop("`n_pcs` must be at least 1.", call. = FALSE)
  }

  score_cols <- sprintf("PC%02d", seq_len(n_pcs))
  scores_dt <- data.table::as.data.table(pca$x[, seq_len(n_pcs), drop = FALSE])
  data.table::setnames(scores_dt, score_cols)

  list(
    pca_model = pca,
    embedding_cols = embedding_cols,
    scores = data.table::data.table(filtered_dt[, .(image_name, id, label)], scores_dt),
    score_cols = score_cols,
    keep_rows = keep_rows
  )
}

cluster_reduced_embeddings <- function(
  reduced_tbl,
  feature_cols = NULL,
  k_min = 2,
  k_max = 10,
  nstart = 20,
  seed = 1
) {
  required_packages <- c("data.table", "cluster")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      sprintf("Install required packages before running: %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }

  dt <- data.table::as.data.table(reduced_tbl)
  if (is.null(feature_cols)) {
    feature_cols <- grep("^PC[0-9]+$", names(dt), value = TRUE)
  }

  if (length(feature_cols) == 0) {
    stop("No reduced-dimension feature columns found. Pass `feature_cols` explicitly.", call. = FALSE)
  }

  feature_dt <- dt[, ..feature_cols]
  feature_mat <- as.matrix(feature_dt)
  storage.mode(feature_mat) <- "double"

  keep_rows <- stats::complete.cases(feature_mat)
  if (!any(keep_rows)) {
    stop("No complete reduced-dimension rows available for clustering.", call. = FALSE)
  }

  filtered_dt <- dt[keep_rows]
  filtered_mat <- feature_mat[keep_rows, , drop = FALSE]
  n_obs <- nrow(filtered_mat)

  if (n_obs < 2) {
    stop("At least two rows are required for clustering.", call. = FALSE)
  }

  k_min <- max(2L, as.integer(k_min))
  k_max <- min(as.integer(k_max), n_obs - 1L)
  if (k_max < k_min) {
    stop("No valid cluster counts available. Check `k_min`, `k_max`, and the number of rows.", call. = FALSE)
  }

  dist_mat <- stats::dist(filtered_mat)
  candidate_ks <- seq.int(k_min, k_max)

  candidate_metrics <- lapply(candidate_ks, function(k) {
    set.seed(seed)
    km <- stats::kmeans(filtered_mat, centers = k, nstart = nstart)
    sil <- cluster::silhouette(km$cluster, dist_mat)
    list(
      k = k,
      model = km,
      mean_silhouette = mean(sil[, "sil_width"])
    )
  })

  metrics_dt <- data.table::rbindlist(lapply(candidate_metrics, function(x) {
    data.table::data.table(
      k = x$k,
      mean_silhouette = x$mean_silhouette
    )
  }))

  best_idx <- which.max(metrics_dt$mean_silhouette)
  best_result <- candidate_metrics[[best_idx]]
  centers <- best_result$model$centers[best_result$model$cluster, , drop = FALSE]
  cluster_distance <- sqrt(rowSums((filtered_mat - centers) ^ 2))

  clustered_dt <- data.table::copy(filtered_dt)
  clustered_dt[, cluster := factor(best_result$model$cluster)]
  clustered_dt[, cluster_distance := cluster_distance]

  list(
    clustered_tbl = clustered_dt,
    metrics = metrics_dt,
    best_k = best_result$k,
    best_mean_silhouette = best_result$mean_silhouette,
    kmeans_model = best_result$model,
    feature_cols = feature_cols,
    keep_rows = keep_rows
  )
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

read_tiff_array <- function(path) {
  tiff::readTIFF(path, native = FALSE, all = FALSE)
}

read_mask_tiff <- function(path) {
  tiff::readTIFF(path, as.is = TRUE)
}

normalize_channel <- function(x, contrast_quantile) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  lo <- min(x, na.rm = TRUE)
  hi <- stats::quantile(x, probs = contrast_quantile, na.rm = TRUE, names = FALSE, type = 7)
  if (!is.finite(hi) || hi <= lo) {
    hi <- max(x, na.rm = TRUE)
  }
  if (!is.finite(hi) || hi <= lo) {
    return(array(0, dim = dim(x)))
  }
  out <- (x - lo) / (hi - lo)
  out[out < 0] <- 0
  out[out > 1] <- 1
  out
}

to_rgb_image <- function(image, contrast_quantile) {
  dims <- dim(image)

  if (is.null(dims) || length(dims) == 2) {
    channel <- normalize_channel(image, contrast_quantile)
    return(array(rep(channel, 3), dim = c(nrow(channel), ncol(channel), 3)))
  }

  if (length(dims) != 3) {
    stop("Unsupported image dimensionality.", call. = FALSE)
  }

  if (dims[[3]] == 1) {
    channel <- normalize_channel(image[, , 1], contrast_quantile)
    return(array(rep(channel, 3), dim = c(dims[[1]], dims[[2]], 3)))
  }

  channels <- min(dims[[3]], 3)
  out <- array(0, dim = c(dims[[1]], dims[[2]], 3))
  for (idx in seq_len(channels)) {
    out[, , idx] <- normalize_channel(image[, , idx], contrast_quantile)
  }
  if (channels == 2) {
    out[, , 3] <- out[, , 2]
  }
  out
}

mask_boundaries <- function(mask) {
  storage.mode(mask) <- "integer"
  nr <- nrow(mask)
  nc <- ncol(mask)
  boundary <- matrix(FALSE, nrow = nr, ncol = nc)

  if (nr > 1) {
    diff_down <- mask[-1, , drop = FALSE] != mask[-nr, , drop = FALSE]
    boundary[-1, ] <- boundary[-1, ] | diff_down
    boundary[-nr, ] <- boundary[-nr, ] | diff_down
  }
  if (nc > 1) {
    diff_right <- mask[, -1, drop = FALSE] != mask[, -nc, drop = FALSE]
    boundary[, -1] <- boundary[, -1] | diff_right
    boundary[, -nc] <- boundary[, -nc] | diff_right
  }

  boundary & mask > 0L
}

blend_boundary <- function(image_rgb, boundary, color_rgb, alpha) {
  out <- image_rgb
  for (idx in seq_len(3)) {
    channel <- out[, , idx]
    channel[boundary] <- (1 - alpha) * channel[boundary] + alpha * color_rgb[[idx]]
    out[, , idx] <- channel
  }
  out
}

write_rgb_tiff <- function(image_rgb, path) {
  tiff::writeTIFF(image_rgb, path, bits.per.sample = 8L, compression = "LZW")
}

default_cluster_palette <- function(n) {
  grDevices::hcl.colors(n, palette = "Dark 3")
}

render_cluster_overlays <- function(
  clustered_tbl,
  segmentation_manifest = "manifests/dev_subset_segmentation_images.csv",
  out_dir = "data/dev_subset_cluster_overlays",
  cluster_col = "cluster",
  image_col = "image_name",
  label_col = "label",
  boundary_alpha = 0.85,
  contrast_quantile = 0.995,
  palette = NULL
) {
  required_packages <- c("data.table", "tiff")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      sprintf("Install required packages before running: %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }

  cluster_dt <- data.table::as.data.table(clustered_tbl)
  required_cols <- c(image_col, label_col, cluster_col)
  missing_cols <- setdiff(required_cols, names(cluster_dt))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("Missing required cluster columns: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  project_root <- resolve_project_root()
  segmentation_manifest <- project_path(project_root, segmentation_manifest)
  manifest_dt <- data.table::fread(segmentation_manifest)
  manifest_dt[, raw_relpath := project_path(project_root, normalize_relpath(raw_relpath))]
  manifest_dt[, mask_relpath := project_path(project_root, normalize_relpath(mask_relpath))]
  manifest_dt[, image_name := filename]

  cluster_dt <- cluster_dt[, .(
    image_name = get(image_col),
    label = as.integer(get(label_col)),
    cluster = as.character(get(cluster_col))
  )]

  cluster_dt <- unique(cluster_dt[!is.na(image_name) & !is.na(label) & !is.na(cluster)])
  if (nrow(cluster_dt) == 0) {
    stop("No valid clustered rows available for overlay rendering.", call. = FALSE)
  }

  render_dt <- merge(
    cluster_dt,
    manifest_dt[, .(image_name, raw_relpath, mask_relpath)],
    by = "image_name",
    all.x = TRUE
  )

  missing_images <- unique(render_dt[is.na(raw_relpath) | is.na(mask_relpath), image_name])
  if (length(missing_images) > 0) {
    stop(
      sprintf(
        "No segmentation manifest entries found for: %s",
        paste(missing_images, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  cluster_levels <- sort(unique(render_dt$cluster))
  if (is.null(palette)) {
    palette <- default_cluster_palette(length(cluster_levels))
  }
  if (length(palette) < length(cluster_levels)) {
    stop("`palette` must provide at least one color per cluster.", call. = FALSE)
  }

  cluster_colors <- setNames(lapply(seq_along(cluster_levels), function(i) {
    as.numeric(grDevices::col2rgb(palette[[i]])[, 1] / 255)
  }), cluster_levels)

  ensure_dir(out_dir)
  written_paths <- vector("list", length = 0L)

  for (current_image_name in unique(render_dt$image_name)) {
    image_rows <- render_dt[image_name == current_image_name]
    image <- read_tiff_array(image_rows$raw_relpath[[1]])
    mask <- read_mask_tiff(image_rows$mask_relpath[[1]])
    if (length(dim(mask)) != 2) {
      stop(sprintf("Expected 2D label mask: %s", image_rows$mask_relpath[[1]]), call. = FALSE)
    }

    image_rgb <- to_rgb_image(image, contrast_quantile)
    base_boundary <- mask_boundaries(mask)

    for (i in seq_len(nrow(image_rows))) {
      object_boundary <- base_boundary & mask == image_rows$label[[i]]
      if (!any(object_boundary)) {
        next
      }
      image_rgb <- blend_boundary(
        image_rgb,
        object_boundary,
        cluster_colors[[image_rows$cluster[[i]]]],
        boundary_alpha
      )
    }

    stem <- sub("\\.[^.]+$", "", basename(current_image_name))
    out_path <- file.path(out_dir, sprintf("%s_cluster_overlay.tif", stem))
    write_rgb_tiff(image_rgb, out_path)
    written_paths[[length(written_paths) + 1L]] <- data.table::data.table(
      image_name = current_image_name,
      out_path = normalize_relpath(out_path)
    )
  }

  list(
    overlays = data.table::rbindlist(written_paths, use.names = TRUE, fill = TRUE),
    cluster_colors = setNames(palette[seq_along(cluster_levels)], cluster_levels)
  )
}
