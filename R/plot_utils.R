# Convert a parent-child edge table into a phylo object.
# Inputs: edge data.frame and optional column names / edge-length callback.
# Output: `ape::phylo` tree.
edge_df_to_phylo <- function(df, id_col = "passage_id", parent_col = "passage_from", edge_length = NULL) {
  library(dplyr)
  library(ape)

  df <- df %>%
    dplyr::transmute(id = .data[[id_col]], parent = .data[[parent_col]]) %>%
    dplyr::distinct()

  stopifnot(!anyDuplicated(df$id))
  stopifnot(!any(is.na(df$id)))

  roots <- df$id[is.na(df$parent) | !(df$parent %in% df$id)]
  if (length(roots) != 1) stop("Expected exactly one root; found: ", paste(roots, collapse = ", "))

  child_map <- split(df$id[!is.na(df$parent)], df$parent[!is.na(df$parent)])
  child_map <- lapply(child_map, unique)

  quote_label <- function(x) {
    paste0("'", gsub("'", "", x, fixed = TRUE), "'")
  }

  build_newick <- function(node) {
    kids <- child_map[[node]]
    if (is.null(kids) || !length(kids)) {
      return(quote_label(node))
    }

    kids <- sort(kids)
    paste0(
      "(",
      paste(vapply(kids, build_newick, character(1)), collapse = ","),
      ")",
      quote_label(node)
    )
  }

  phy <- ape::read.tree(text = paste0(build_newick(roots[1]), ";"))
  phy$tip.label <- gsub("^'|'$", "", phy$tip.label)
  phy$tip.label <- gsub(" ", "_", phy$tip.label, fixed = TRUE)
  if (!is.null(phy$node.label)) {
    phy$node.label <- gsub("^'|'$", "", phy$node.label)
    phy$node.label <- gsub(" ", "_", phy$node.label, fixed = TRUE)
  }
  if (!is.null(edge_length)) {
    el <- edge_length(df)
    if (length(el) == nrow(phy$edge)) {
      phy$edge.length <- as.numeric(el)
    }
  }
  phy
}

# Build a rooted layout directly from the passage graph.
# Input: edge table with `passage_id` and `passage_from`.
# Output: list with node coordinates, edge coordinates, and root id.
layout_passage_tree <- function(edge_df) {
  edge_df <- unique(edge_df[, c("passage_id", "passage_from"), drop = FALSE])
  edge_df$passage_id <- as.character(edge_df$passage_id)
  edge_df$passage_from <- as.character(edge_df$passage_from)

  if (anyDuplicated(edge_df$passage_id)) {
    stop("layout_passage_tree() requires unique passage_id values.")
  }
  if (any(is.na(edge_df$passage_id))) {
    stop("layout_passage_tree() requires non-missing passage_id values.")
  }

  node_ids <- unique(c(edge_df$passage_id, edge_df$passage_from))
  node_ids <- node_ids[!is.na(node_ids)]
  roots <- unique(edge_df$passage_id[is.na(edge_df$passage_from) | !(edge_df$passage_from %in% edge_df$passage_id)])
  if (length(roots) != 1L) {
    stop("Expected exactly one root; found: ", paste(roots, collapse = ", "))
  }

  child_map <- split(edge_df$passage_id[!is.na(edge_df$passage_from)], edge_df$passage_from[!is.na(edge_df$passage_from)])
  child_map <- lapply(child_map, function(children) sort(unique(as.character(children))))

  depths <- stats::setNames(rep(NA_real_, length(node_ids)), node_ids)
  assign_depths <- function(node_id, depth_value) {
    depths[[node_id]] <<- depth_value
    children <- child_map[[node_id]]
    if (is.null(children) || !length(children)) {
      return(invisible(NULL))
    }
    for (child_id in children) {
      assign_depths(child_id, depth_value + 1)
    }
  }
  assign_depths(roots[1], 0)

  leaf_order <- 0
  y_pos <- stats::setNames(rep(NA_real_, length(node_ids)), node_ids)
  assign_y <- function(node_id) {
    children <- child_map[[node_id]]
    if (is.null(children) || !length(children)) {
      leaf_order <<- leaf_order + 1
      y_pos[[node_id]] <<- leaf_order
      return(y_pos[[node_id]])
    }

    child_y <- vapply(children, assign_y, numeric(1))
    y_pos[[node_id]] <<- mean(child_y)
    y_pos[[node_id]]
  }
  assign_y(roots[1])

  node_df <- data.frame(
    label = node_ids,
    x = as.numeric(depths[node_ids]),
    y = as.numeric(y_pos[node_ids]),
    stringsAsFactors = FALSE
  )

  edge_plot_df <- edge_df[!is.na(edge_df$passage_from), , drop = FALSE]
  if (nrow(edge_plot_df)) {
    edge_plot_df$x_parent <- node_df$x[match(edge_plot_df$passage_from, node_df$label)]
    edge_plot_df$y_parent <- node_df$y[match(edge_plot_df$passage_from, node_df$label)]
    edge_plot_df$x_child <- node_df$x[match(edge_plot_df$passage_id, node_df$label)]
    edge_plot_df$y_child <- node_df$y[match(edge_plot_df$passage_id, node_df$label)]
  }

  list(
    nodes = node_df,
    edges = edge_plot_df,
    root_id = roots[1]
  )
}

# Plot one connected passage tree with nodes colored by image sample count.
# Inputs: `tree` list from the workflow and optional `component_name`.
# Output: `ggplot` plot object.
plot_passage_tree <- function(tree, component_name = NULL) {
  library(ggplot2)

  tree_layout <- layout_passage_tree(tree$tree_edges)
  node_df <- tree_layout$nodes

  sample_counts <- table(as.character(tree$tree_samples$passage_id))
  node_df$image_count <- as.integer(sample_counts[node_df$label])
  node_df$image_count[is.na(node_df$image_count)] <- 0L
  node_df$image_count_bin <- cut(
    node_df$image_count,
    breaks = c(-Inf, 0, 1, 5, Inf),
    labels = c("0", "1", "2-5", "6+"),
    right = TRUE
  )
  node_df$image_count_bin <- factor(node_df$image_count_bin, levels = c("0", "1", "2-5", "6+"))
  condition_df <- tree$tree_conditions
  if (is.null(condition_df)) {
    condition_df <- data.frame(
      passage_id = character(),
      condition = character(),
      stringsAsFactors = FALSE
    )
  }
  node_df$condition <- condition_df$condition[match(node_df$label, condition_df$passage_id)]
  node_df$condition[is.na(node_df$condition) | !nzchar(node_df$condition)] <- "Unannotated"
  node_df$stressor_class <- condition_df$stressor_class[match(node_df$label, condition_df$passage_id)]
  node_df$stressor_class[is.na(node_df$stressor_class) | !nzchar(node_df$stressor_class)] <- "none"

  parent_segments <- do.call(
    rbind,
    lapply(split(tree_layout$edges, tree_layout$edges$passage_from), function(parent_edges) {
      data.frame(
        x = parent_edges$x_parent[1],
        xend = parent_edges$x_parent[1],
        y = min(parent_edges$y_child),
        yend = max(parent_edges$y_child),
        stringsAsFactors = FALSE
      )
    })
  )
  if (is.null(parent_segments)) {
    parent_segments <- data.frame(x = numeric(), xend = numeric(), y = numeric(), yend = numeric())
  }
  tree_layout$edges$condition <- node_df$condition[match(tree_layout$edges$passage_id, node_df$label)]
  tree_layout$edges$condition[is.na(tree_layout$edges$condition) | !nzchar(tree_layout$edges$condition)] <- "Unannotated"

  ggplot() +
    geom_segment(
      data = parent_segments,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.35,
      color = "grey70",
      lineend = "round"
    ) +
    geom_segment(
      data = tree_layout$edges,
      aes(x = x_parent, xend = x_child, y = y_child, yend = y_child, color = condition),
      linewidth = 1.5,
      lineend = "round"
    ) +
    geom_point(
      data = node_df,
      aes(x = x, y = y, fill = image_count_bin),
      size = 2.4,
      alpha = 0.95,
      shape = 21,
      stroke = 0.25,
      color = "grey20"
    ) +
    geom_point(
      data = node_df[node_df$stressor_class %in% "other", , drop = FALSE],
      aes(x = x, y = y),
      shape = 8,
      size = 2.8,
      stroke = 0.6,
      color = "black"
    ) +
    geom_point(
      data = node_df[node_df$stressor_class %in% "reversine", , drop = FALSE],
      aes(x = x, y = y),
      shape = 4,
      size = 3,
      stroke = 0.9,
      color = "red3"
    ) +
    scale_color_brewer("media\ncondition", palette = "Set2", na.translate = FALSE) +
    scale_fill_manual(
      "image\nsamples",
      values = c("0" = "#440154FF", "1" = "#31688EFF", "2-5" = "#35B779FF", "6+" = "#FDE725FF"),
      drop = FALSE
    ) +
    scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
    theme_void() +
    ggtitle(if (is.null(component_name)) "Connected passage tree" else paste("Connected passage tree:", component_name))
}

# Plot a legacy experimental lineage colored by fitted growth rate.
# Input: `lineage_df_subset` edge table with metadata columns used in the plot.
# Output: `ggtree` plot object.
experiment_lineage_plot <- function(lineage_df_subset){
  library(dplyr)
  library(ggplot2)
  library(ggtree)

  phy <- edge_df_to_phylo(lineage_df_subset)

  meta <- lineage_df_subset %>%
    dplyr::transmute(
      label = as.character(passage_id),
      g = ifelse(is.finite(g), g, NA_real_),
      karyotyped = as.logical(karyotyped),
      has_flow = as.logical(has_flow)
    )

  p <- ggtree(phy, layout="circular", size=0.3)
  p$data <- dplyr::left_join(p$data, meta, by="label")

  p <- p +
    geom_point2(
      data = ~ dplyr::filter(.x, karyotyped),
      aes(shape = "karyotyped"), color = "red", size = 3.2
    ) +
    geom_point2(
      data = ~ dplyr::filter(.x, has_flow),
      aes(shape = "has_flow"), color = "green", size = 3.2
    ) +
    geom_point2(
      aes(color = g),
      shape = 16, size = 1.4
    ) +
    scale_color_viridis_c("growth\nrate", na.value = "grey70", trans = "log") +
    theme_void() +
    scale_shape_discrete("") +
    ggtitle(paste("Hypoxia experiment data:", unique(lineage_df_subset$label_value)))
  return(p)
}
