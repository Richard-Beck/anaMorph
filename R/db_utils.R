# Read database credentials from a simple KEY=value file.
# Input: `filepath` path to a text file containing HOST, DBNAME, USER, PASSWORD.
# Output: named character vector of credential values.
load_db_vars <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("Error: File '%s' does not exist.", filepath))
  }
  
  lines <- tryCatch(readLines(filepath, warn = FALSE), 
                    error = function(e) stop("Error: Failed to read file '", filepath, "'. ", e$message))
  
  if (length(lines) == 0) {
    stop(sprintf("Error: File '%s' is empty.", filepath))
  }
  
  if (any(!grepl("^[A-Za-z_][A-Za-z0-9_]*=.+$", lines))) {
    stop("Error: All lines must be in the format KEY=value, with valid variable names.")
  }
  
  kv_pairs <- strsplit(lines, "=", fixed = TRUE)
  keys <- vapply(kv_pairs, `[`, "", 1)
  values <- vapply(kv_pairs, `[`, "", 2)
  vars <- setNames(values, keys)
  
  required_keys <- c("HOST", "DBNAME", "USER", "PASSWORD")
  missing_keys <- setdiff(required_keys, names(vars))
  if (length(missing_keys) > 0) {
    stop("Error: Missing required keys: ", paste(missing_keys, collapse = ", "))
  }
  
  return(vars)
}

# Fetch per-sample karyotype vectors from the database.
# Input: `cloneIds` character vector of sample ids to query.
# Output: tibble with columns `id` and list-column `karyotype`.
get_karyotyping <- function(cloneIds){
  require(DBI)
  require(RMariaDB)
  dbvars <- load_db_vars("db_creds.txt")
  db <- dbConnect(
    MariaDB(),
    host = dbvars["HOST"],
    user = dbvars["USER"],
    password = dbvars["PASSWORD"],
    dbname = dbvars["DBNAME"]
  )
  
  origin_clause <- paste0("'", cloneIds, "'", collapse = ", ")
  q <- paste0("SELECT * FROM Perspective WHERE origin IN (", origin_clause, ") AND whichPerspective='GenomePerspective'")
  
  rs <- dbSendQuery(db, q)
  res <- dbFetch(rs)
  dbClearResult(rs) 
  dbDisconnect(db)
  
  kvecs <- lapply(res$profile, function(p) {
    raw_vec <- p
    readBin(raw_vec, what = "double", n = length(raw_vec) / 8, endian = "big")
  })
  
  df <- tibble::tibble(
    id = res$origin,
    karyotype = kvecs
  )
  df <- df[!res$hasChildren==1,] ## this SHOULD exclude the population average entry (...?)
  df
}

# Collapse raw passaging records to a passage-level ancestry graph.
# Input: `x` data.frame containing `id`, `event`, `passaged_from_id1`, and `media`.
# Output: list with `passaging` edges and `samples` mapping raw ids to collapsed passages.
collapse_db <- function(x) {
  required_cols <- c("id", "event", "passaged_from_id1", "media")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop("collapse_db() missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  if (anyDuplicated(x$id)) {
    stop("collapse_db() requires unique values in x$id.")
  }

  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x$id <- as.character(x$id)
  x$event <- as.character(x$event)
  x$passaged_from_id1 <- as.character(x$passaged_from_id1)
  x$media <- as.character(x$media)
  x$passaged_from_id1[is.na(x$passaged_from_id1) | x$passaged_from_id1 == ""] <- NA_character_
  x$media[is.na(x$media) | x$media == ""] <- NA_character_

  valid_events <- c("seeding", "harvest")
  bad_event_rows <- !(x$event %in% valid_events) | is.na(x$event)
  if (any(bad_event_rows)) {
    warning(
      "collapse_db() stripped ", sum(bad_event_rows),
      " rows with unsupported event values."
    )
    x <- x[!bad_event_rows, , drop = FALSE]
  }

  missing_parent_rows <- !is.na(x$passaged_from_id1) & !(x$passaged_from_id1 %in% x$id)
  if (any(missing_parent_rows)) {
    warning(
      "collapse_db() stripped ", sum(missing_parent_rows),
      " rows whose parent id was absent from x$id."
    )
    x <- x[!missing_parent_rows, , drop = FALSE]
  }

  root_rows <- is.na(x$passaged_from_id1)
  is_passage_row <- x$event == "seeding" | root_rows
  valid_passage_ids <- x$id[is_passage_row]
  parent_lookup <- setNames(x$passaged_from_id1, x$id)

  resolve_passage_id <- function(id, parent_map, valid_ids) {
    cur <- id
    visited <- character()
    while (!is.null(cur) && !is.na(cur)) {
      if (cur %in% visited) {
        stop("collapse_db() detected a cycle while resolving passage ids.")
      }
      if (cur %in% valid_ids) return(cur)
      visited <- c(visited, cur)
      cur <- unname(parent_map[cur])
    }
    NA_character_
  }

  sample_passage_id <- vapply(
    x$id,
    resolve_passage_id,
    character(1),
    parent_map = parent_lookup,
    valid_ids = valid_passage_ids
  )

  unresolved_rows <- is.na(sample_passage_id)
  if (any(unresolved_rows)) {
    warning(
      "collapse_db() stripped ", sum(unresolved_rows),
      " rows that could not be resolved to any upstream passage_id."
    )
    x <- x[!unresolved_rows, , drop = FALSE]
    sample_passage_id <- sample_passage_id[!unresolved_rows]
    root_rows <- root_rows[!unresolved_rows]
    is_passage_row <- is_passage_row[!unresolved_rows]
    valid_passage_ids <- valid_passage_ids[valid_passage_ids %in% x$id]
    parent_lookup <- setNames(x$passaged_from_id1, x$id)
  }

  samples <- data.frame(
    passage_id = sample_passage_id,
    samples = x$id,
    media_id = x$media,
    stringsAsFactors = FALSE
  )

  passaging_rows <- x[is_passage_row, , drop = FALSE]
  parent_passage_id <- vapply(
    passaging_rows$passaged_from_id1,
    resolve_passage_id,
    character(1),
    parent_map = parent_lookup,
    valid_ids = valid_passage_ids
  )
  parent_passage_id[is.na(passaging_rows$passaged_from_id1)] <- NA_character_

  passaging <- data.frame(
    passage_id = passaging_rows$id,
    passage_from = unname(parent_passage_id),
    stringsAsFactors = FALSE
  )

  list(
    passaging = passaging,
    samples = samples
  )
}

# Subset the media table to the media ids used in a sample set.
# Inputs: `samples` vector or data.frame containing media ids, `media_tbl` raw media table.
# Output: reduced media data.frame with `media_id` replacing `id`.
build_media_col <- function(samples, media_tbl) {
  if (is.data.frame(samples)) {
    if (!("media_id" %in% names(samples))) {
      stop("build_media_col() requires a 'media_id' column when samples is a data.frame.")
    }
    media_ids <- unique(as.character(samples$media_id))
  } else {
    media_ids <- unique(as.character(samples))
  }

  media_ids <- media_ids[!is.na(media_ids) & nzchar(media_ids)]

  if (!is.data.frame(media_tbl) || !("id" %in% names(media_tbl))) {
    stop("build_media_col() requires media_tbl to be a data.frame containing an 'id' column.")
  }

  out <- media_tbl[media_tbl$id %in% media_ids, , drop = FALSE]
  if (!nrow(out)) {
    out <- media_tbl[FALSE, , drop = FALSE]
  }

  names(out)[names(out) == "id"] <- "media_id"
  keep_cols <- vapply(out, function(col) length(unique(col)) > 1L, logical(1))
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out[, keep_cols, drop = FALSE]
}

# Extract ancestor-descendant intervals between karyotyped passages in one tree.
# Input: `tree` list containing sampled passages, tree edges, and tree conditions.
# Output: list of interval records with start/end passages and path conditions.
extract_karyotyped_intervals <- function(tree) {
  required_tree_fields <- c("sampled_passage_ids", "tree_edges", "tree_conditions")
  missing_tree_fields <- setdiff(required_tree_fields, names(tree))
  if (length(missing_tree_fields) > 0) {
    stop(
      "extract_karyotyped_intervals() missing required tree fields: ",
      paste(missing_tree_fields, collapse = ", ")
    )
  }

  edge_df <- tree$tree_edges
  if (!all(c("passage_id", "passage_from") %in% names(edge_df))) {
    stop("extract_karyotyped_intervals() requires tree$tree_edges to contain 'passage_id' and 'passage_from'.")
  }

  condition_df <- tree$tree_conditions
  if (!all(c("passage_id", "condition") %in% names(condition_df))) {
    stop("extract_karyotyped_intervals() requires tree$tree_conditions to contain 'passage_id' and 'condition'.")
  }

  parent_of <- setNames(as.character(edge_df$passage_from), as.character(edge_df$passage_id))
  condition_of <- setNames(as.character(condition_df$condition), as.character(condition_df$passage_id))
  karyotyped_ids <- unique(as.character(tree$sampled_passage_ids))

  find_direct_karyotyped_ancestor <- function(node_id, parent_map, sampled_ids) {
    cur <- unname(parent_map[node_id])
    while (!is.na(cur) && !(cur %in% sampled_ids)) {
      cur <- unname(parent_map[cur])
    }
    if (is.na(cur)) {
      return(NULL)
    }
    cur
  }

  get_path_between <- function(start_id, end_id, parent_map) {
    out <- end_id
    cur <- end_id
    visited <- character()

    while (!identical(cur, start_id)) {
      if (cur %in% visited) {
        stop("extract_karyotyped_intervals() detected a cycle while traversing ancestry.")
      }
      visited <- c(visited, cur)
      cur <- unname(parent_map[cur])
      if (is.null(cur) || is.na(cur)) {
        stop("Passage ", start_id, " is not an ancestor of ", end_id, ".")
      }
      out <- c(cur, out)
    }

    out
  }

  is_valid_direct_interval <- function(start_id, end_id, parent_map, sampled_ids) {
    path_ids <- tryCatch(
      get_path_between(start_id, end_id, parent_map),
      error = function(e) NULL
    )
    if (is.null(path_ids)) {
      return(FALSE)
    }

    intermediate_ids <- setdiff(path_ids, c(start_id, end_id))
    !any(intermediate_ids %in% sampled_ids)
  }

  build_interval <- function(start_id, end_id, start_alias = NULL) {
    if (is.null(start_alias)) {
      path_ids <- get_path_between(start_id, end_id, parent_of)
    } else {
      path_ids <- c(start_id, get_path_between(start_alias, end_id, parent_of))
    }
    path_conditions <- condition_of[path_ids]
    path_conditions[is.na(path_conditions) | !nzchar(path_conditions)] <- "Unannotated"

    list(
      start = start_id,
      end = end_id,
      start_alias = start_alias,
      n_collapsed_passages = if (is.null(start_alias)) 0L else 1L,
      conditions = table(path_conditions)
    )
  }

  direct_intervals <- lapply(karyotyped_ids, function(end_id) {
    direct_start <- find_direct_karyotyped_ancestor(end_id, parent_of, karyotyped_ids)
    if (!is.null(direct_start)) {
      return(build_interval(direct_start, end_id, start_alias = NULL))
    }
    NULL
  })
  direct_intervals <- direct_intervals[!vapply(direct_intervals, is.null, logical(1))]

  valid_direct_starts <- unique(vapply(direct_intervals, `[[`, character(1), "start"))
  endpoints_with_direct_start <- unique(vapply(direct_intervals, `[[`, character(1), "end"))

  alias_intervals <- list()
  alias_idx <- 1L
  proxy_starts <- setdiff(karyotyped_ids, valid_direct_starts)
  for (proxy_start in proxy_starts) {
    start_alias <- unname(parent_of[proxy_start])
    if (is.null(start_alias) || is.na(start_alias)) {
      next
    }

    candidate_ends <- setdiff(karyotyped_ids, c(proxy_start, endpoints_with_direct_start))
    for (end_id in candidate_ends) {
      if (!is_valid_direct_interval(start_alias, end_id, parent_of, karyotyped_ids)) {
        next
      }

      alias_intervals[[alias_idx]] <- build_interval(
        start_id = proxy_start,
        end_id = end_id,
        start_alias = start_alias
      )
      alias_idx <- alias_idx + 1L
    }
  }

  intervals <- c(direct_intervals, alias_intervals)
  if (!length(intervals)) {
    return(intervals)
  }

  interval_keys <- vapply(intervals, function(interval) {
    paste(interval$start, interval$end, if (is.null(interval$start_alias)) "NULL" else interval$start_alias, sep = "||")
  }, character(1))
  intervals[!duplicated(interval_keys)]
}

# Collapse a path-level condition table to one interval label.
# Input: `condition_counts` named table of conditions along an interval.
# Output: scalar character condition label.
resolve_interval_condition <- function(condition_counts) {
  condition_names <- names(condition_counts)
  if (!length(condition_names)) {
    return("unknown")
  }

  if (length(condition_names) == 1L) {
    return(condition_names[1])
  }

  if (length(condition_names) == 2L && "control" %in% condition_names) {
    return(setdiff(condition_names, "control")[1])
  }

  "unknown"
}

# Add karyotype matrices and metadata to tree-derived intervals.
# Inputs: `connected_trees` list and optional `karyotypes` table with `id` and `karyotype`.
# Output: flat list of simulation-ready interval records.
build_simulation_intervals <- function(connected_trees, karyotypes = NULL) {
  if (is.null(names(connected_trees))) {
    names(connected_trees) <- as.character(seq_along(connected_trees))
  }

  empty_karyotype <- matrix(NA_real_, nrow = 1, ncol = 1)
  keep_chr_idx <- 1:22

  karyotype_matrix_for_ids <- function(sample_ids, karyotype_tbl) {
    sample_ids <- unique(as.character(sample_ids))
    sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]

    if (!length(sample_ids) || is.null(karyotype_tbl)) {
      return(empty_karyotype)
    }

    if (!all(c("id", "karyotype") %in% names(karyotype_tbl))) {
      stop("build_simulation_intervals() requires karyotypes to contain 'id' and 'karyotype' columns.")
    }

    rows <- karyotype_tbl[as.character(karyotype_tbl$id) %in% sample_ids, , drop = FALSE]
    if (!nrow(rows)) {
      return(empty_karyotype)
    }

    kmat <- do.call(rbind, lapply(rows$karyotype, function(k) as.numeric(unlist(k))))
    kmat <- as.matrix(kmat)
    if (ncol(kmat) < max(keep_chr_idx)) {
      stop("build_simulation_intervals() expected at least 22 chromosome columns in karyotypes.")
    }
    kmat <- kmat[, keep_chr_idx, drop = FALSE]
    rownames(kmat) <- make.unique(as.character(rows$id))
    kmat
  }

  intervals_by_tree <- lapply(seq_along(connected_trees), function(i) {
    tree_name <- names(connected_trees)[i]
    tree <- connected_trees[[i]]
    tree_intervals <- extract_karyotyped_intervals(tree)

    lapply(tree_intervals, function(interval) {
      start_ids <- tree$tree_samples$samples[tree$tree_samples$passage_id == interval$start]
      end_ids <- tree$tree_samples$samples[tree$tree_samples$passage_id == interval$end]

      interval$connected_tree <- tree_name
      interval$condition <- resolve_interval_condition(interval$conditions)
      interval$elapsed_passages <- as.integer(sum(interval$conditions))
      interval$start_ids <- unique(as.character(start_ids))
      interval$end_ids <- unique(as.character(end_ids))
      interval$start_karyotype <- karyotype_matrix_for_ids(interval$start_ids, karyotypes)
      interval$end_karyotype <- karyotype_matrix_for_ids(interval$end_ids, karyotypes)
      interval
    })
  })

  unlist(intervals_by_tree, recursive = FALSE)
}

# Group intervals into shared parameter sets.
# Inputs: `simulation_intervals` list and `divisions_per_passage` scalar integer.
# Output: named list of interval groups keyed by connected tree and condition.
group_simulation_intervals <- function(simulation_intervals, divisions_per_passage = 5L) {
  if (!length(simulation_intervals)) {
    return(list())
  }

  intervals <- lapply(simulation_intervals, function(interval) {
    interval$n_divisions <- as.integer(interval$elapsed_passages * divisions_per_passage)
    interval$group_id <- paste(interval$connected_tree, interval$condition, sep = "||")
    interval
  })

  split(intervals, vapply(intervals, `[[`, character(1), "group_id"))
}

# Recover one raw lineage between two samples and annotate simple growth fits.
# Inputs: `row` containing `first`/`last` ids, `data` raw passaging table.
# Output: data.frame of lineage rows with fitted growth columns added.
recover_lineage <- function(row, data) {
  current_id <- row$last
  lineage <- c(current_id)
  
  while (!is.na(current_id) && current_id != row$first) {
    parent <- data$passaged_from_id1[which(data$id == current_id)]
    if (length(parent) == 0 || is.na(parent)) {
      stop(paste("Could not trace from", row$last, "to", row$first))
    }
    lineage <- c(parent, lineage)
    current_id <- parent
  }
  
  if (current_id != row$first) {
    stop(paste("Could not find", row$first, "starting from", row$last))
  }
  
  seeding_lineage <- lineage[data$event[match(lineage, data$id)] == "seeding"]
  
  filtered_x <- do.call(rbind, lapply(seq_along(seeding_lineage), function(i) {
    dfi <- data[data$id %in% seeding_lineage[i] | data$passaged_from_id1 %in% seeding_lineage[i], ]
    dfi$adjPass <- stringr::str_pad(i, width = 2)
    dfi$date <- as.Date(dfi$date)
    dfi$num_date <- as.numeric(as.Date(dfi$date))
    dfi$num_date <- dfi$num_date - min(dfi$num_date)
    dfi$intercept <- NaN
    dfi$g <- NaN
    if (nrow(dfi) < 2) return(dfi)
    if(sum(!is.na(dfi$correctedCount))<2) return(dfi)
    fit <- lm(log(pmax(1, dfi$correctedCount)) ~ dfi$num_date)
    dfi$intercept <- exp(coef(fit)[1])
    dfi$g <- coef(fit)[2]
    dfi
  }))
  
  filtered_x$label_value    <- row$label
  filtered_x$sublabel_value <- row$label2
  return(filtered_x)
}


# Fill missing lineage segments for a SUM-159 subset in legacy exploratory analysis.
# Inputs: ploidy substring, current lineage data.frame, and raw passaging table `x`.
# Output: filtered data.frame of lineage rows to append.
fill_lineage_gaps <- function(ploidy_substr,df,x){
  g <- x[x$cellLine=="SUM-159",c("id","passaged_from_id1")]
  
  
  ## perhaps split 2N 4N.
  karyotyped_lineages <- lapply(cloneIds[grepl(ploidy_substr,cloneIds)],function(id){
    ids <- c()
    while(!is.na(id)){
      ids <- c(id,ids)
      id <- g$passaged_from_id1[g$id==id]
    }
    ids
  })
  
  init_ids <- sapply(filters,'[[',"first")
  init_ids <- init_ids[grepl(ploidy_substr,init_ids)]
  
  imaged_lineages <- lapply(init_ids,function(id){
    ids <- c()
    while(!is.na(id)){
      ids <- c(id,ids)
      id <- g$passaged_from_id1[g$id==id]
    }
    ids
  })
  
  all_lineages <- c(imaged_lineages,karyotyped_lineages)
  n_lineages <- length(all_lineages)
  
  uids <- table(unlist(all_lineages))
  n_remove <- length(uids[uids==n_lineages])-1
  
  all_lineages <- lapply(all_lineages, function(i){
    i[-c(1:n_remove)]
  })
  
  x <- x[x$id%in%unlist(all_lineages),]
  x <- x[!x$id%in%df$id,]
  x$label_value <- gsub("_","",ploidy_substr)
  x
}
