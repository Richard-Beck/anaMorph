parse_args <- function(args) {
  parsed <- list(
    segmentation_manifest = "manifests/dev_subset_segmentation_images.csv",
    out_dir = "data/dev_subset_overlays",
    boundary_color = "#FF3B30",
    boundary_alpha = "0.85",
    contrast_quantile = "0.995"
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

  required <- c("segmentation_manifest", "out_dir")
  missing <- required[!nzchar(vapply(required, function(x) {
    val <- parsed[[x]]
    if (is.null(val)) "" else as.character(val)
  }, character(1)))]
  if (length(missing) > 0) {
    stop(sprintf("Missing required args: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  parsed$boundary_alpha <- as.numeric(parsed$boundary_alpha)
  parsed$contrast_quantile <- as.numeric(parsed$contrast_quantile)
  parsed
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

normalize_relpath <- function(path) {
  gsub("\\\\", "/", path)
}

require_tiff <- function() {
  if (!requireNamespace("tiff", quietly = TRUE)) {
    stop(
      "Package 'tiff' is required. Install it with install.packages('tiff').",
      call. = FALSE
    )
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

hex_to_rgb <- function(hex_color) {
  rgb_vals <- grDevices::col2rgb(hex_color) / 255
  as.numeric(rgb_vals[, 1])
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

main <- function() {
  require_tiff()
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  manifest <- read.csv(args$segmentation_manifest, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(manifest) == 0) {
    stop("No rows found in segmentation manifest.", call. = FALSE)
  }

  ensure_dir(args$out_dir)
  overlay_color <- hex_to_rgb(args$boundary_color)

  written <- 0L
  for (row_idx in seq_len(nrow(manifest))) {
    row <- manifest[row_idx, ]
    image <- read_tiff_array(row$raw_relpath[[1]])
    mask <- read_mask_tiff(row$mask_relpath[[1]])
    if (length(dim(mask)) != 2) {
      stop(sprintf("Expected 2D label mask: %s", row$mask_relpath[[1]]), call. = FALSE)
    }

    image_rgb <- to_rgb_image(image, args$contrast_quantile)
    boundary <- mask_boundaries(mask)
    overlay <- blend_boundary(image_rgb, boundary, overlay_color, args$boundary_alpha)

    stem <- sub("\\.[^.]+$", "", basename(row$filename[[1]]))
    out_path <- file.path(args$out_dir, sprintf("%s_overlay.tif", stem))
    write_rgb_tiff(overlay, out_path)
    written <- written + 1L
  }

  cat(sprintf("Wrote %d overlays to %s\n", written, normalize_relpath(args$out_dir)))
}

main()
