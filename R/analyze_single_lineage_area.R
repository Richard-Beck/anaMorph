parse_args <- function(args) {
  parsed <- list(
    lineage_rds = "core_data/lineages.Rds",
    object_sizes_rds = "data/all_images/object_sizes.rds",
    out_dir = "data/lineage_area"
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

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

parse_mixed_datetime <- function(x) {
  x <- as.character(x)
  out <- as.POSIXct(rep(NA_character_, length(x)), tz = "UTC")
  formats <- c(
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%dT%H:%M:%SZ",
    "%m/%d/%Y %H:%M",
    "%m/%d/%Y %H:%M:%S"
  )

  for (fmt in formats) {
    missing <- is.na(out) & !is.na(x) & nzchar(x)
    if (!any(missing)) next
    parsed <- as.POSIXct(x[missing], format = fmt, tz = "UTC")
    out[missing] <- parsed
  }

  out
}

safe_quantile <- function(x, prob) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7))
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  stats::weighted.mean(x[ok], w[ok])
}

weighted_var_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  mu <- stats::weighted.mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

weighted_skew_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 3) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  mu <- stats::weighted.mean(x, w)
  v <- sum(w * (x - mu)^2) / sum(w)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  m3 <- sum(w * (x - mu)^3) / sum(w)
  m3 / (v^(3 / 2))
}

fit_two_normal_mixture <- function(x, max_iter = 200, tol = 1e-6) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 50) {
    mu <- stats::median(x, na.rm = TRUE)
    sigma <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(sigma) || sigma <= 0) sigma <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sigma) || sigma <= 0) sigma <- 1
    return(list(
      pi = c(0.05, 0.95),
      mu = c(mu - sigma, mu),
      sigma = c(sigma, sigma),
      posterior_small = rep(0, length(x))
    ))
  }

  q <- stats::quantile(x, probs = c(0.15, 0.65), na.rm = TRUE, names = FALSE)
  s0 <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s0) || s0 <= 0) s0 <- 1

  pi_k <- c(0.15, 0.85)
  mu_k <- c(q[[1]], q[[2]])
  sigma_k <- c(s0 * 0.5, s0)
  loglik_prev <- -Inf

  for (iter in seq_len(max_iter)) {
    dens1 <- pi_k[[1]] * stats::dnorm(x, mean = mu_k[[1]], sd = sigma_k[[1]])
    dens2 <- pi_k[[2]] * stats::dnorm(x, mean = mu_k[[2]], sd = sigma_k[[2]])
    denom <- pmax(dens1 + dens2, .Machine$double.eps)
    gamma1 <- dens1 / denom
    gamma2 <- dens2 / denom

    n1 <- sum(gamma1)
    n2 <- sum(gamma2)
    pi_k <- c(n1, n2) / length(x)
    mu_k <- c(sum(gamma1 * x) / n1, sum(gamma2 * x) / n2)
    sigma_k <- c(
      sqrt(sum(gamma1 * (x - mu_k[[1]])^2) / n1),
      sqrt(sum(gamma2 * (x - mu_k[[2]])^2) / n2)
    )
    sigma_k[!is.finite(sigma_k) | sigma_k < 1e-3] <- 1e-3

    loglik <- sum(log(denom))
    if (abs(loglik - loglik_prev) < tol) break
    loglik_prev <- loglik
  }

  small_idx <- which.min(mu_k)
  dens_small <- pi_k[[small_idx]] * stats::dnorm(x, mean = mu_k[[small_idx]], sd = sigma_k[[small_idx]])
  dens_other <- pi_k[[3 - small_idx]] * stats::dnorm(x, mean = mu_k[[3 - small_idx]], sd = sigma_k[[3 - small_idx]])
  posterior_small <- dens_small / pmax(dens_small + dens_other, .Machine$double.eps)

  list(
    pi = pi_k,
    mu = mu_k,
    sigma = sigma_k,
    posterior_small = posterior_small,
    small_component = small_idx
  )
}

make_formula <- function(response, data) {
  terms <- character()

  unique_time <- length(unique(stats::na.omit(data$time_since_passage_start_hours)))
  if (unique_time >= 4) {
    df_time <- min(3, unique_time - 1)
    terms <- c(terms, sprintf("splines::ns(time_since_passage_start_hours, df = %d)", df_time))
  } else if (unique_time >= 2) {
    terms <- c(terms, "time_since_passage_start_hours")
  }

  unique_conf <- length(unique(stats::na.omit(data$log_confluence)))
  if (unique_conf >= 4) {
    df_conf <- min(3, unique_conf - 1)
    terms <- c(terms, sprintf("splines::ns(log_confluence, df = %d)", df_conf))
  } else if (unique_conf >= 2) {
    terms <- c(terms, "log_confluence")
  }

  if ("field" %in% names(data) && length(unique(stats::na.omit(data$field))) >= 2) {
    terms <- c(terms, "field")
  }

  if (!length(terms)) {
    stats::as.formula(sprintf("%s ~ 1", response))
  } else {
    stats::as.formula(sprintf("%s ~ %s", response, paste(terms, collapse = " + ")))
  }
}

extract_lineage <- function(lineages, lineage_id) {
  if (lineage_id %in% names(lineages)) {
    return(lineages[[lineage_id]])
  }

  hit <- vapply(lineages, function(x) identical(as.character(x$lineage_id), as.character(lineage_id)), logical(1))
  if (!any(hit)) {
    stop(sprintf("Lineage id not found: %s", lineage_id), call. = FALSE)
  }

  lineages[[which(hit)[[1]]]]
}

build_lineage_cell_table <- function(lineage_obj, object_sizes) {
  passage_positions <- stats::setNames(seq_along(lineage_obj$passage_ids), lineage_obj$passage_ids)
  image_meta_list <- lapply(names(lineage_obj$passages), function(passage_id) {
    passage_dt <- data.table::as.data.table(lineage_obj$passages[[passage_id]])
    if (!nrow(passage_dt)) return(NULL)

    keep_cols <- intersect(
      c("id", "filename", "field", "date", "areaOccupied_um2", "cellCount", "condition_label", "event"),
      names(passage_dt)
    )
    image_meta <- unique(passage_dt[, ..keep_cols], by = "filename")
    image_meta[, passage_id := as.character(passage_id)]
    image_meta[, lineage_id := as.character(lineage_obj$lineage_id)]
    image_meta[, passage_position := unname(passage_positions[passage_id])]
    image_meta[, image_key := sub("\\.[^.]+$", "", basename(filename))]
    image_meta
  })

  image_meta <- data.table::rbindlist(image_meta_list, fill = TRUE)
  if (!nrow(image_meta)) {
    stop("No manifest rows found in lineage object.", call. = FALSE)
  }

  image_meta[, image_time := parse_mixed_datetime(date)]
  image_meta[, passage_start_time := suppressWarnings(min(image_time, na.rm = TRUE)), by = passage_id]
  image_meta[!is.finite(as.numeric(passage_start_time)), passage_start_time := as.POSIXct(NA, tz = "UTC")]
  image_meta[, time_since_passage_start_hours := as.numeric(difftime(image_time, passage_start_time, units = "hours"))]
  image_meta[, confluence_um2 := suppressWarnings(as.numeric(areaOccupied_um2))]
  image_meta[, log_confluence := log1p(pmax(confluence_um2, 0))]

  matched_keys <- intersect(image_meta$image_key, names(object_sizes))
  if (!length(matched_keys)) {
    stop("No image keys from this lineage matched names(object_sizes).", call. = FALSE)
  }

  image_meta <- image_meta[image_key %in% matched_keys]
  cell_tables <- lapply(seq_len(nrow(image_meta)), function(i) {
    row <- image_meta[i]
    object_dt <- data.table::as.data.table(object_sizes[[row$image_key]])
    if (!nrow(object_dt)) return(NULL)

    object_dt[, label := suppressWarnings(as.integer(label))]
    object_dt[, area_px := suppressWarnings(as.numeric(area_px))]
    object_dt <- object_dt[is.finite(area_px) & area_px > 0]
    if (!nrow(object_dt)) return(NULL)

    object_dt[, log_area := log(area_px)]
    object_dt[, lineage_id := row$lineage_id]
    object_dt[, passage_id := row$passage_id]
    object_dt[, passage_position := row$passage_position]
    object_dt[, image_id := row$id]
    object_dt[, image_key := row$image_key]
    object_dt[, filename := row$filename]
    object_dt[, field := if ("field" %in% names(row)) row$field else NA_character_]
    object_dt[, image_time := row$image_time]
    object_dt[, time_since_passage_start_hours := row$time_since_passage_start_hours]
    object_dt[, confluence_um2 := row$confluence_um2]
    object_dt[, log_confluence := row$log_confluence]
    object_dt[, cell_count := if ("cellCount" %in% names(row)) suppressWarnings(as.numeric(row$cellCount)) else NA_real_]
    object_dt[, condition_label := if ("condition_label" %in% names(row)) row$condition_label else NA_character_]
    object_dt
  })

  data.table::rbindlist(cell_tables, fill = TRUE)
}

summarize_group <- function(dt, group_cols) {
  dt[, .(
    n_objects = .N,
    effective_n = sum(p_real, na.rm = TRUE),
    nuisance_fraction = mean(p_small, na.rm = TRUE),
    median_z = stats::median(z, na.rm = TRUE),
    iqr_z = stats::IQR(z, na.rm = TRUE),
    sd_z = stats::sd(z, na.rm = TRUE),
    skew_z = weighted_skew_safe(z, p_real),
    q95_z = safe_quantile(z, 0.95),
    q99_z = safe_quantile(z, 0.99),
    frac_z_gt_2 = mean(z > 2, na.rm = TRUE),
    frac_z_gt_3 = mean(z > 3, na.rm = TRUE),
    mean_time_since_passage_start_hours = mean(time_since_passage_start_hours, na.rm = TRUE),
    mean_log_confluence = mean(log_confluence, na.rm = TRUE)
  ), by = group_cols]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- resolve_project_root()
lineage_rds <- project_path(project_root, args$lineage_rds)
object_sizes_rds <- project_path(project_root, args$object_sizes_rds)
out_dir <- project_path(project_root, file.path(args$out_dir, paste0("lineage_", args$lineage_id)))

ensure_dir(out_dir)

lineages <- readRDS(lineage_rds)
object_sizes <- readRDS(object_sizes_rds)
lineage_obj <- extract_lineage(lineages, args$lineage_id)
cell_dt <- build_lineage_cell_table(lineage_obj, object_sizes)

if (!nrow(cell_dt)) {
  stop("No object rows available for this lineage.", call. = FALSE)
}

mixture_fit <- fit_two_normal_mixture(cell_dt$log_area)
cell_dt[, p_small := mixture_fit$posterior_small]
cell_dt[, p_real := pmax(1 - p_small, 1e-4)]

model_dt <- cell_dt[is.finite(log_area) & is.finite(time_since_passage_start_hours) & is.finite(log_confluence)]
if (nrow(model_dt) < 50) {
  stop("Too few modeled cells after covariate filtering.", call. = FALSE)
}

location_formula <- make_formula("log_area", model_dt)
location_fit <- MASS::rlm(location_formula, data = model_dt, weights = model_dt$p_real, maxit = 100)
cell_dt[, mu_hat := as.numeric(stats::predict(location_fit, newdata = cell_dt))]

model_dt[, mu_hat := as.numeric(stats::predict(location_fit, newdata = model_dt))]
model_dt[, abs_resid := pmax(abs(log_area - mu_hat), 1e-4)]
scale_formula <- make_formula("log(abs_resid)", model_dt)
scale_fit <- stats::lm(scale_formula, data = model_dt, weights = model_dt$p_real)
cell_dt[, sigma_hat := exp(as.numeric(stats::predict(scale_fit, newdata = cell_dt)))]

global_sigma <- stats::mad(model_dt$log_area - model_dt$mu_hat, constant = 1.4826, na.rm = TRUE)
if (!is.finite(global_sigma) || global_sigma <= 0) {
  global_sigma <- stats::sd(model_dt$log_area - model_dt$mu_hat, na.rm = TRUE)
}
if (!is.finite(global_sigma) || global_sigma <= 0) {
  global_sigma <- 1
}

cell_dt[!is.finite(sigma_hat) | sigma_hat <= 0, sigma_hat := global_sigma]
cell_dt[, z := (log_area - mu_hat) / sigma_hat]

image_summary <- summarize_group(
  cell_dt,
  c("lineage_id", "passage_id", "passage_position", "image_id", "image_key", "filename", "field")
)

passage_summary <- summarize_group(
  cell_dt,
  c("lineage_id", "passage_id", "passage_position")
)

passage_summary <- merge(
  passage_summary,
  unique(cell_dt[, .(
    passage_id,
    first_image_time = if (all(is.na(image_time))) as.POSIXct(NA, tz = "UTC") else min(image_time, na.rm = TRUE),
    mean_confluence_um2 = mean(confluence_um2, na.rm = TRUE)
  ), by = passage_id]),
  by = "passage_id",
  all.x = TRUE,
  sort = FALSE
)

out_rds <- file.path(out_dir, "lineage_area_analysis.rds")
out_cells_csv <- file.path(out_dir, "cell_residuals.csv")
out_images_csv <- file.path(out_dir, "image_summary.csv")
out_passages_csv <- file.path(out_dir, "passage_summary.csv")

saveRDS(list(
  lineage_id = lineage_obj$lineage_id,
  passage_ids = lineage_obj$passage_ids,
  mixture_fit = mixture_fit,
  location_formula = deparse(location_formula),
  scale_formula = deparse(scale_formula),
  location_coefficients = stats::coef(location_fit),
  scale_coefficients = stats::coef(scale_fit),
  cells = cell_dt,
  image_summary = image_summary,
  passage_summary = passage_summary
), out_rds)

utils::write.csv(cell_dt, out_cells_csv, row.names = FALSE, na = "")
utils::write.csv(image_summary, out_images_csv, row.names = FALSE, na = "")
utils::write.csv(passage_summary, out_passages_csv, row.names = FALSE, na = "")

cat(sprintf("Lineage: %s\n", lineage_obj$lineage_id))
cat(sprintf("Passages: %d\n", length(lineage_obj$passage_ids)))
cat(sprintf("Images matched to object sizes: %d\n", length(unique(cell_dt$image_key))))
cat(sprintf("Objects analyzed: %d\n", nrow(cell_dt)))
cat(sprintf("Wrote analysis RDS: %s\n", out_rds))
