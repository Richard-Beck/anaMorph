# Build the ancestor path from the component root to a passage.
# Inputs: `passage_id` scalar character, `parent_map` named character vector of parent ids.
# Output: character vector ordered from root to `passage_id`.
find_passage_ancestors <- function(passage_id, parent_map) {
  out <- character()
  cur <- passage_id
  while (!is.null(cur) && !is.na(cur)) {
    out <- c(cur, out)
    if (!(cur %in% names(parent_map))) {
      break
    }
    cur <- unname(parent_map[cur])
  }
  out
}

# Find the shared prefix across several root-to-node paths.
# Inputs: `paths` list of character vectors.
# Output: character vector giving the shared prefix.
common_path_prefix <- function(paths) {
  if (!length(paths)) {
    return(character())
  }
  
  min_len <- min(vapply(paths, length, integer(1)))
  prefix <- character()
  for (i in seq_len(min_len)) {
    vals <- vapply(paths, `[`, character(1), i)
    if (length(unique(vals)) != 1L) {
      break
    }
    prefix <- c(prefix, vals[1])
  }
  prefix
}

# Drop the upstream part of a path so it starts at a requested node.
# Inputs: `path` character vector, `node_id` scalar character that must appear in `path`.
# Output: suffix of `path` beginning at `node_id`.
trim_path_to_node <- function(path, node_id) {
  node_idx <- match(node_id, path)
  if (is.na(node_idx)) {
    stop("Node ", node_id, " was not present in the supplied path.")
  }
  path[node_idx:length(path)]
}

# Build a connected passage tree trimmed to the MRCA of sampled passages.
# Inputs: component passage ids, precomputed path map, collapsed passaging table, and sample map.
# Output: list describing one connected tree plus sampled ids and trimmed paths.
build_component_tree <- function(component_passage_ids, path_map, passaging_tbl, sample_map) {
  component_paths <- path_map[component_passage_ids]
  common_prefix <- common_path_prefix(component_paths)
  if (!length(common_prefix)) {
    stop("Could not determine a shared ancestry path for component.")
  }
  
  mrca_id <- tail(common_prefix, 1)
  trimmed_paths <- lapply(component_paths, trim_path_to_node, node_id = mrca_id)
  subtree_nodes <- sort(unique(unlist(trimmed_paths)))
  subtree_edges <- passaging_tbl[
    passaging_tbl$passage_id %in% subtree_nodes,
    c("passage_id", "passage_from"),
    drop = FALSE
  ]
  subtree_edges$passage_from[!(subtree_edges$passage_from %in% subtree_nodes)] <- NA_character_
  subtree_samples <- sample_map[sample_map$passage_id %in% subtree_nodes, , drop = FALSE]
  
  list(
    component_root_id = common_prefix[1],
    mrca_id = mrca_id,
    mrca_depth = length(common_prefix),
    sampled_passage_ids = sort(component_passage_ids),
    sampled_ids = sort(unique(subtree_samples$samples[subtree_samples$passage_id %in% component_passage_ids])),
    tree_passage_ids = subtree_nodes,
    tree_edges = subtree_edges[order(subtree_edges$passage_id), , drop = FALSE],
    tree_samples = subtree_samples[order(subtree_samples$passage_id, subtree_samples$samples), , drop = FALSE],
    trimmed_paths = trimmed_paths
  )
}

# Convert sampled ids into connected trees over the collapsed passage graph.
# Inputs: `db_col` from `collapse_db()`, `sample_ids` character vector of karyotyped ids.
# Output: named list with `sample_map`, `connected_trees`, and `component_summary`.
build_connected_trees <- function(db_col, sample_ids) {
  sample_ids <- unique(as.character(sample_ids))
  sample_ids <- intersect(sample_ids, db_col$samples$samples)
  
  sample_map <- db_col$samples[match(sample_ids, db_col$samples$samples), , drop = FALSE]
  sample_map <- unique(sample_map[, c("samples", "passage_id"), drop = FALSE])
  sample_map <- sample_map[order(sample_map$passage_id, sample_map$samples), , drop = FALSE]
  
  passaging_tbl <- db_col$passaging
  parent_of <- setNames(passaging_tbl$passage_from, passaging_tbl$passage_id)
  unmapped_passage_ids <- setdiff(unique(sample_map$passage_id), passaging_tbl$passage_id)
  if (length(unmapped_passage_ids) > 0) {
    warning(
      length(unmapped_passage_ids),
      " sampled passage_id values are absent from db_col$passaging and will be skipped."
    )
  }
  
  sampled_passage_ids <- setdiff(unique(sample_map$passage_id), unmapped_passage_ids)
  passage_paths <- setNames(
    lapply(sampled_passage_ids, find_passage_ancestors, parent_map = parent_of),
    sampled_passage_ids
  )
  root_id <- vapply(passage_paths, function(path) path[1], character(1))
  connected_passage_sets <- split(sampled_passage_ids, root_id)
  
  connected_trees <- lapply(
    connected_passage_sets,
    build_component_tree,
    path_map = passage_paths,
    passaging_tbl = passaging_tbl,
    sample_map = sample_map
  )
  
  component_summary <- do.call(rbind, lapply(seq_along(connected_trees), function(i) {
    tree <- connected_trees[[i]]
    data.frame(
      component = names(connected_trees)[i],
      component_root_id = tree$component_root_id,
      mrca_id = tree$mrca_id,
      trimmed_pre_mrca = tree$component_root_id != tree$mrca_id,
      mrca_depth = tree$mrca_depth,
      n_sampled_passages = length(tree$sampled_passage_ids),
      n_sampled_ids = length(tree$sampled_ids),
      n_tree_passages = length(tree$tree_passage_ids),
      stringsAsFactors = FALSE
    )
  }))
  
  list(
    sample_map = sample_map,
    connected_trees = connected_trees,
    component_summary = component_summary
  )
}


# Resolve one representative media id per passage using the most common sample annotation.
# Inputs: `passage_ids` character vector, `samples_tbl` with `passage_id` and `media_id`.
# Output: data.frame with one row per passage and columns `passage_id` and `media_id`.
resolve_passage_media <- function(passage_ids, samples_tbl) {
  media_ids <- lapply(passage_ids, function(id) {
    samples_tbl$media_id[samples_tbl$passage_id %in% id]
  })
  resolved_media <- vapply(media_ids, function(ids) {
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (!length(ids)) {
      return(NA_character_)
    }
    agg <- table(ids)
    names(agg)[which.max(agg)]
  }, character(1))
  
  data.frame(
    passage_id = passage_ids,
    media_id = unname(resolved_media),
    stringsAsFactors = FALSE
  )
}

# Resolve passage-level stressor status from any associated sample media.
# Inputs: `passage_ids`, `samples_tbl` with `passage_id`/`media_id`, and raw `media_tbl`.
# Output: data.frame with `passage_id`, `has_stressor`, and `stressor_class`.
resolve_passage_stressor <- function(passage_ids, samples_tbl, media_tbl) {
  media_lookup <- setNames(as.character(media_tbl$Stressor), as.character(media_tbl$id))
  media_lookup[is.na(media_lookup)] <- ""

  stressor_info <- lapply(passage_ids, function(id) {
    media_ids <- as.character(samples_tbl$media_id[samples_tbl$passage_id %in% id])
    stressors <- unname(media_lookup[media_ids])
    stressors <- unique(stressors[!is.na(stressors) & nzchar(stressors)])

    if (!length(stressors)) {
      return(c(has_stressor = "FALSE", stressor_class = "none"))
    }
    if ("Reversine" %in% stressors) {
      return(c(has_stressor = "TRUE", stressor_class = "reversine"))
    }
    c(has_stressor = "TRUE", stressor_class = "other")
  })

  data.frame(
    passage_id = passage_ids,
    has_stressor = vapply(stressor_info, function(x) identical(x[["has_stressor"]], "TRUE"), logical(1)),
    stressor_class = vapply(stressor_info, `[[`, character(1), "stressor_class"),
    stringsAsFactors = FALSE
  )
}

# Label media rows with coarse experimental conditions using hand-written rules.
# Inputs: `media_tbl` from `build_media_col()`.
# Output: `media_tbl` with an added `condition` column.
annotate_media_conditions <- function(media_tbl) {
  filters <- list(
    !is.na(media_tbl$EnergySource2) & media_tbl$EnergySource2_pct < 100 & media_tbl$EnergySource2=="Phosphates",
    media_tbl$oxygen_pct < 20.5,
    media_tbl$EnergySource == "L_glutamine" &
      media_tbl$EnergySource_nM < 2000000 &
      !is.na(media_tbl$EnergySource_nM),
    media_tbl$EnergySource == "Glucose" &
      media_tbl$EnergySource_nM < 1110150 &
      !is.na(media_tbl$EnergySource_nM)
  )
  filters <- do.call(cbind, filters)
  filters[is.na(filters)] <- FALSE
  
  media_tbl$condition <- apply(filters, 1, function(row_flags) {
    if (sum(row_flags) == 0L) {
      return("control")
    }
    deprivations <- c("phosphate", "oxygen", "glutamine", "glucose")[row_flags]
    paste(deprivations, collapse = "_")
  })
  
  media_tbl
}

# Attach per-passage condition annotations to each connected tree.
# Inputs: `connected_trees`, `db_samples` from `db_col$samples`, and raw `media_tbl`.
# Output: list with updated `connected_trees`, `passage_media`, and annotated `media_col`.
assign_conditions_to_trees <- function(connected_trees, db_samples, media_tbl) {
  all_ids <- unique(unlist(lapply(connected_trees, `[[`, "tree_passage_ids")))
  passage_media <- resolve_passage_media(all_ids, db_samples)
  passage_stressor <- resolve_passage_stressor(all_ids, db_samples, media_tbl)
  media_col <- build_media_col(unique(passage_media$media_id), media_tbl)
  media_col <- annotate_media_conditions(media_col)
  passage_media$condition <- media_col$condition[match(passage_media$media_id, media_col$media_id)]
  passage_media$has_stressor <- passage_stressor$has_stressor[match(passage_media$passage_id, passage_stressor$passage_id)]
  passage_media$stressor_class <- passage_stressor$stressor_class[match(passage_media$passage_id, passage_stressor$passage_id)]
  
  connected_trees <- lapply(connected_trees, function(tree) {
    tree$tree_conditions <- passage_media[
      match(tree$tree_passage_ids, passage_media$passage_id),
      ,
      drop = FALSE
    ]
    tree$tree_conditions <- tree$tree_conditions[order(tree$tree_conditions$passage_id), , drop = FALSE]
    tree
  })
  
  list(
    connected_trees = connected_trees,
    passage_media = passage_media,
    media_col = media_col
  )
}
