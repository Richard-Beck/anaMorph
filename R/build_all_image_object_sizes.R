parse_args <- function(args) {
  parsed <- list(
    embedding_dir = "data/all_images/embeddings",
    out_rds = "data/all_images/object_sizes.rds"
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
        dir.exists(file.path(candidate, "data"))) {
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

ensure_parent_dir <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- resolve_project_root()
embedding_dir <- project_path(project_root, args$embedding_dir)
out_rds <- project_path(project_root, args$out_rds)

if (!dir.exists(embedding_dir)) {
  stop(sprintf("Embedding directory not found: %s", embedding_dir), call. = FALSE)
}

embedding_files <- list.files(
  embedding_dir,
  pattern = "_embeddings\\.csv$",
  full.names = TRUE
)

if (length(embedding_files) == 0) {
  stop(sprintf("No embedding CSVs found in %s", embedding_dir), call. = FALSE)
}

object_sizes <- pbapply::pblapply(embedding_files, function(path) {
  data.table::fread(path, select = c("label", "area_px"))
})

names(object_sizes) <- sub("_embeddings\\.csv$", "", basename(embedding_files))

ensure_parent_dir(out_rds)
saveRDS(object_sizes, out_rds)

test_read <- readRDS(out_rds)
cat(sprintf("Read %d embedding CSVs from %s\n", length(object_sizes), embedding_dir))
cat(sprintf("Wrote object sizes RDS: %s\n", out_rds))
cat(sprintf("Verified RDS contains %d list entries\n", length(test_read)))
