source("results/baps_palette.R")
#!/usr/bin/env Rscript
# Circular P. aeruginosa phylogeny adapted directly from
# /Users/hc/ST58_project-main/st131/plotcase_phylogeny/code.R.
# Only the input tree and metadata fields are changed for this project.

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggtree)
  library(ggplot2)
  library(ggnewscale)
  library(cowplot)
  library(readr)
  library(dplyr)
})

set.seed(131)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)][1])
if (!nzchar(script) || identical(script, "NA")) script <- file.path(getwd(), "analysis_3_1", "plot_pseudomonas_circular_tree.R")
work_dir <- dirname(normalizePath(script, mustWork = FALSE))
root <- normalizePath(file.path(work_dir, "..", ".."), mustWork = TRUE)
out_dir <- file.path(work_dir, "fan_tree_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Input files: Verticall recombination-filtered tree and the 3.1 annotation table.
tree_file <- file.path(root, "tree", "verticall-rm-recombination-tree", "verticall.newick")
metadata_file <- file.path(work_dir, "sample_annotation_3_1.tsv")
if (!file.exists(metadata_file)) metadata_file <- file.path(root, "results", "3.1_analysis", "sample_annotation_3_1.tsv")
stopifnot(file.exists(tree_file), file.exists(metadata_file))

clean_annotation <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(z) | tolower(z) %in% c("", "missing", "unknown", "unknow", "na", "n/a")] <- NA_character_
  z
}
collapse_rare <- function(x, threshold = 10L) {
  z <- clean_annotation(x)
  counts <- table(z, useNA = "no")
  z[!is.na(z) & z %in% names(counts)[counts < threshold]] <- "Others"
  z
}
ordered_categories <- function(x) {
  z <- sort(table(x, useNA = "no"), decreasing = TRUE)
  c(setdiff(names(z), "Others"), intersect("Others", names(z)))
}

metadata <- readr::read_tsv(metadata_file, show_col_types = FALSE,
                            na = c("", "NA", "<NA>"))
tree_input <- ape::read.tree(tree_file)
if (anyDuplicated(tree_input$tip.label) || anyDuplicated(metadata$sample_id)) {
  stop("Tree tip labels and metadata sample_id values must be unique")
}
if (is.null(tree_input$edge.length) || any(!is.finite(tree_input$edge.length))) {
  stop("Verticall tree must contain finite branch lengths")
}
negative_edges <- sum(tree_input$edge.length < 0)
if (negative_edges > 0L) {
  warning("Verticall tree contains ", negative_edges,
          " negative branch lengths; they are retained as supplied for the distance audit.")
}
if (!setequal(tree_input$tip.label, metadata$sample_id)) {
  stop("Tree and metadata sample sets do not match")
}

# Same midpoint-rooting procedure as code.R. Midpoint rooting changes display
# orientation only; pairwise tip distances are audited below.
tree <- phangorn::midpoint(tree_input)
stopifnot(is.rooted(tree), setequal(tree$tip.label, tree_input$tip.label))
distance_before <- cophenetic(tree_input)
distance_after <- cophenetic(tree)[rownames(distance_before), colnames(distance_before)]
distance_error <- max(abs(distance_after - distance_before), na.rm = TRUE)
diameter_idx <- arrayInd(which.max(distance_before), dim(distance_before))[1, ]
diameter_tips <- rownames(distance_before)[diameter_idx]
root_depth <- setNames(node.depth.edgelength(tree)[seq_len(Ntip(tree))], tree$tip.label)
midpoint_error <- abs(diff(root_depth[diameter_tips]))
write.tree(tree, file.path(out_dir, "Pseudomonas.verticall.midpoint_rooted.nwk"))

# Match metadata by the tree-driven tip order, exactly as in code.R.
meta <- metadata[match(tree$tip.label, metadata$sample_id), , drop = FALSE]
rownames(meta) <- NULL
ann <- data.frame(Name = tree$tip.label, check.names = FALSE)
ann$fastBAPS <- clean_annotation(meta$fastBAPS_cluster)
ann$ST <- collapse_rare(meta$st)
ann$Country <- clean_annotation(meta$country)
ann$Year <- clean_annotation(meta$year_group)
ann$Serotype <- clean_annotation(meta$serotype)
ann$`Source proxy` <- clean_annotation(meta$source_proxy)

# Four outer rings show whether each isolate carries a strict LoF event in the
# focal genes.  The LoF rings are binary and use gene-specific colours.
lof_file <- file.path(root, "results", "3.6_LoF_analysis", "results", "LoF_events_with_metadata.tsv")
if (!file.exists(lof_file)) lof_file <- file.path(root, "3.6_LoF_analysis", "results", "LoF_events_with_metadata.tsv")
lof <- readr::read_tsv(lof_file, show_col_types = FALSE)
lof_genes <- c("oprD", "mucA", "mexB", "lasR")
lof_present <- setNames(lapply(lof_genes, function(g) {
  unique(as.character(lof$sample[lof$gene_key == g]))
}), lof_genes)
for (g in lof_genes) ann[[paste0("LoF_", g)]] <- ifelse(ann$Name %in% lof_present[[g]], g, "Absent")
ring_names <- c("ST", "Country", "Year", "Serotype", paste0("LoF: ", lof_genes))
write.csv(ann, file.path(out_dir, "Pseudomonas.plot_annotations.csv"), row.names = FALSE, na = "")

# Explicit palettes are kept separate for each ring, as in code.R.
categorical_palette <- c("#8DD3C7", "#E8CE54", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#52A69B",
  "#B79C62", "#E777AA", "#5273B0", "#B76B4C", "#638C50", "#B1B8C2",
  "#665A9B", "#36B6C5", "#CE9DCC", "#A3C94A")
palettes <- lapply(ring_names, function(nm) {
  if (startsWith(nm, "LoF: ")) {
    gene <- sub("^LoF: ", "", nm)
    return(c(Absent = "white", setNames(c("#D55E00", "#0072B2", "#009E73", "#CC79A7")[match(gene, lof_genes)], gene)))
  }
  categories <- ordered_categories(ann[[nm]])
  cols <- if (length(categories) <= length(categorical_palette)) categorical_palette[seq_along(categories)] else
    grDevices::hcl.colors(length(categories), palette = "Dynamic")
  p <- setNames(cols, categories)
  if ("Others" %in% categories) p["Others"] <- "#A1A1A1"
  p
})
names(palettes) <- ring_names
fast_levels <- ordered_categories(ann$fastBAPS)
fast_palette <- setNames(BAPS_COLORS[paste0("BAPS", fast_levels)], fast_levels)

# Fan tree construction copied from code.R. Branch lengths are explicitly
# retained; no equal-length transformation is applied.
rotation <- 230
arc <- 350
p <- ggtree(
  tree,
  layout = "fan",
  open.angle = 360 - arc,
  branch.length = "branch.length",
  linewidth = 0.18,
  color = "#343C40"
) %<+% ann
p <- rotate_tree(p, rotation)
radius <- max(p$data$x, na.rm = TRUE)
tip_data <- p$data[p$data$isTip, , drop = FALSE]
tip_data <- tip_data[match(tree$tip.label, tip_data$label), , drop = FALSE]

# Identify maximal local pure clades for every fastBAPS cluster. A local sector
# is valid when every descendant tip belongs to one cluster; nested pure nodes
# are suppressed so each cluster is represented by the largest non-overlapping
# local sectors available in the supplied topology.
internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
parent_node <- setNames(tree$edge[, 1], tree$edge[, 2])
local_rows <- list()
for (group in fast_levels) {
  group_tips <- ann$Name[!is.na(ann$fastBAPS) & ann$fastBAPS == group]
  pure_by_node <- setNames(vapply(internal_nodes, function(node) {
    descendants <- tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
    length(descendants) >= 2L && all(descendants %in% group_tips)
  }, logical(1)), as.character(internal_nodes))
  pure_nodes <- internal_nodes[pure_by_node[as.character(internal_nodes)]]
  maximal_nodes <- pure_nodes[!vapply(pure_nodes, function(node) {
    parent <- parent_node[as.character(node)]
    !is.na(parent) && as.character(parent) %in% names(pure_by_node) &&
      isTRUE(pure_by_node[[as.character(parent)]])
  }, logical(1))]
  if (length(maximal_nodes)) {
    local_rows[[group]] <- bind_rows(lapply(maximal_nodes, function(node) {
      descendants <- tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
      data.frame(Cluster = group, node = as.integer(node), Samples = length(descendants),
                 Pure = TRUE, stringsAsFactors = FALSE)
    }))
  }
}
local_clades <- bind_rows(local_rows)
if (!nrow(local_clades)) {
  local_clades <- data.frame(Cluster = character(), node = integer(), Samples = integer(), Pure = logical())
}
write.csv(local_clades, file.path(out_dir, "Pseudomonas.fastBAPS_local_clade_sectors.csv"), row.names = FALSE)

# Retain the original whole-cluster MRCA audit as a diagnostic, while using the
# non-overlapping local sectors for the actual fan-tree highlighting.
whole_clade_rows <- lapply(fast_levels, function(group) {
  group_tips <- ann$Name[!is.na(ann$fastBAPS) & ann$fastBAPS == group]
  node <- ape::getMRCA(tree, group_tips)
  descendants <- tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
  data.frame(Cluster = group, node = as.integer(node), Samples = length(group_tips),
             Descendants = length(descendants), Pure = setequal(group_tips, descendants),
             Local_sectors = sum(local_clades$Cluster == group))
})
clades <- bind_rows(whole_clade_rows)
write.csv(clades, file.path(out_dir, "Pseudomonas.fastBAPS_clade_audit.csv"), row.names = FALSE)

if (nrow(local_clades)) {
  branch_layer_count <- length(p$layers)
  p <- p + geom_highlight(data = local_clades, aes(node = node, fill = Cluster),
                          alpha = 0.48, extendto = radius * 1.07) +
    scale_fill_manual(values = BAPS_COLORS, guide = "none", na.translate = FALSE)
  p$layers <- c(p$layers[(branch_layer_count + 1):length(p$layers)], p$layers[seq_len(branch_layer_count)])
  p <- p + ggnewscale::new_scale_fill()
} else {
  p <- p + ggnewscale::new_scale_fill()
}

# Alignment guides and actual tips, following code.R.
ring_width <- radius * 0.052
ring_gap <- radius * 0.006
ring_start <- radius * 1.13
stopifnot(ring_start - ring_width / 2 > radius)
p <- p + geom_segment(
  data = tip_data,
  aes(x = x, xend = ring_start - ring_width / 2, y = y, yend = y),
  inherit.aes = FALSE, color = "#DDE3E5", linewidth = 0.10, linetype = "dotted"
)

ring_centres <- ring_start + (seq_along(ring_names) - 1) * (ring_width + ring_gap)
fast_ring_centre <- ring_start - (ring_width + ring_gap)
fast_ring_values <- ann$fastBAPS[match(tip_data$label, ann$Name)]
p <- p + ggnewscale::new_scale_fill() +
  geom_tile(
    data = data.frame(x = fast_ring_centre, y = tip_data$y, value = fast_ring_values),
    aes(x, y, fill = value), width = ring_width, height = 1,
    color = NA, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = BAPS_COLORS, na.value = "white", na.translate = FALSE,
                    name = "fastBAPS cluster", guide = "none") +
  ggnewscale::new_scale_fill()
for (i in seq_along(ring_names)) {
  nm <- ring_names[i]
  source_nm <- if (startsWith(nm, "LoF: ")) paste0("LoF_", sub("^LoF: ", "", nm)) else nm
  values <- ann[[source_nm]][match(tip_data$label, ann$Name)]
  tile <- data.frame(x = ring_centres[i], y = tip_data$y,
                     value = values, check.names = FALSE)
  p <- p + geom_tile(data = tile, aes(x, y, fill = value),
                     width = ring_width, height = 1, color = NA,
                     inherit.aes = FALSE) +
    scale_fill_manual(values = palettes[[nm]], na.value = "white",
                      na.translate = FALSE, name = nm, guide = "legend")
  if (i < length(ring_names)) p <- p + ggnewscale::new_scale_fill()
}

# FastBAPS is shown at the tips as an explicit point layer; the sectors above
# are drawn only for groups verified as pure clades. Keep the grouping variable
# inside the layer data so every tip, including singleton/non-sector tips, is
# mapped deterministically by ggplot2. A dark outline separates dense tip
# points from branches and highlight sectors at publication-scale resolution.
tip_data$fastBAPS <- ann$fastBAPS[match(tip_data$label, ann$Name)]
p <- p + ggnewscale::new_scale_fill() +
  geom_point(data = tip_data,
    aes(x = x, y = y, fill = fastBAPS), inherit.aes = FALSE,
    shape = 21, color = "#263238", size = 2.0, stroke = 0.35) +
  scale_fill_manual(values = BAPS_COLORS, na.value = "white", na.translate = FALSE,
                    name = "fastBAPS cluster")

# Tangential ring labels in the open fan wedge, copied from code.R.
angular_min <- min(tip_data$y) - 0.5
angular_max <- angular_min + diff(range(tip_data$y)) * 360 / arc
label_y <- max(p$data$y) + (angular_max - max(p$data$y)) * 0.5
label_theta <- (rotation - 90) - 360 * (label_y - angular_min) / (angular_max - angular_min)
label_angle <- (label_theta + 180) %% 360
if (label_angle > 90 && label_angle < 270) label_angle <- label_angle + 180
ring_labels <- data.frame(x = ring_centres, y = label_y, label = ring_names)
p <- p + geom_text(data = ring_labels, aes(x, y, label = label),
  angle = label_angle, size = 2.0, hjust = 0.5, inherit.aes = FALSE, color = "#27343B") +
  scale_y_continuous(limits = c(angular_min, angular_max), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, max(ring_centres) + ring_width * 0.8), expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "right", plot.margin = margin(8, 8, 8, 8),
        plot.background = element_rect(fill = "white", color = NA))

# Validate all ring layers contain one complete cell per tree tip.
built <- ggplot_build(p)
tile_layer_indices <- which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomTile"), logical(1)))
metadata_tile_layer_indices <- tail(tile_layer_indices, length(ring_names))
fast_ring_layer_index <- tile_layer_indices[1]
point_layer_indices <- which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1)))
if (length(point_layer_indices) != 1L) {
  stop("Expected exactly one fastBAPS tip-point layer; found ", length(point_layer_indices))
}
point_data <- built$data[[point_layer_indices[1]]]
rendered_tip_points <- nrow(point_data)
rendered_colored_tip_points <- sum(!is.na(point_data$fill) & point_data$fill != "white")
sector_tip_names <- if (nrow(local_clades)) {
  unique(unlist(lapply(local_clades$node, function(node) {
    tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
  })))
} else character()
nonsector_tip_idx <- !(tip_data$label %in% sector_tip_names)
if (length(tile_layer_indices) != length(ring_names) + 1L) {
  warning("Expected ", length(ring_names) + 1L, " tile layers including fastBAPS strip; found ", length(tile_layer_indices))
}
if (length(metadata_tile_layer_indices) != length(ring_names)) {
  stop("Expected exactly ", length(ring_names), " metadata ring layers")
}
if (nrow(built$data[[fast_ring_layer_index]]) != nrow(tip_data)) {
  stop("FastBAPS tip strip did not render one cell per tip")
}
rendered_fast_strip <- sum(!is.na(built$data[[fast_ring_layer_index]]$fill) &
                           built$data[[fast_ring_layer_index]]$fill != "white")
if (rendered_fast_strip != sum(!is.na(tip_data$fastBAPS))) {
  stop("FastBAPS tip strip did not render all annotated tip fills")
}
rendered_nonsector_colored_tip_points <- sum(nonsector_tip_idx & !is.na(point_data$fill) & point_data$fill != "white")
tip_expected_fill <- unname(fast_palette[as.character(tip_data$fastBAPS)])
if (rendered_tip_points != nrow(tip_data) ||
    rendered_colored_tip_points != sum(!is.na(tip_data$fastBAPS)) ||
    !all(point_data$fill == tip_expected_fill)) {
  stop("FastBAPS tip-point layer did not render all annotated tip fills")
}
ring_cells <- vapply(tile_layer_indices, function(idx) nrow(built$data[[idx]]), integer(1))

png_path <- file.path(out_dir, "Figure_3_1_Pseudomonas_verticall_fan_tree.png")
pdf_path <- file.path(out_dir, "Figure_3_1_Pseudomonas_verticall_fan_tree.pdf")
ggsave(png_path, p, width = 320, height = 320, units = "mm", dpi = 400,
       bg = "white", limitsize = FALSE)
ggsave(pdf_path, p, width = 320, height = 320, units = "mm",
       device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "Figure_3_1_Pseudomonas_verticall_tree_only.png"),
       p, width = 297, height = 297, units = "mm", dpi = 300, bg = "white")

validation <- c(
  sprintf("Tree tips: %s", Ntip(tree)),
  sprintf("Metadata matched: %s", sum(tree$tip.label %in% metadata$sample_id)),
  sprintf("Input negative branch lengths retained: %s", negative_edges),
  sprintf("Midpoint diameter-end depth difference: %.3g", midpoint_error),
  sprintf("Maximum pairwise-distance change after midpoint rooting: %.3g", distance_error),
  sprintf("Rotation: %s degrees; arc: %s degrees; gap: %s degrees", rotation, arc, 360 - arc),
  sprintf("Local fastBAPS clade sectors highlighted: %s across %s clusters", nrow(local_clades), length(unique(local_clades$Cluster))),
  sprintf("Whole-cluster MRCA audit: %s/%s clusters pure", sum(clades$Pure), nrow(clades)),
  sprintf("FastBAPS tip colors mapped: %s/%s", sum(!is.na(tip_data$fastBAPS)), nrow(tip_data)),
  sprintf("Rendered fastBAPS tip points: %s/%s; colored: %s", rendered_tip_points, nrow(tip_data), rendered_colored_tip_points),
  sprintf("Non-sector tips colored explicitly: %s/%s", rendered_nonsector_colored_tip_points, sum(nonsector_tip_idx)),
  sprintf("Metadata ring layers: %s; cells per layer: %s", length(metadata_tile_layer_indices),
          paste(vapply(metadata_tile_layer_indices, function(idx) nrow(built$data[[idx]]), integer(1)), collapse = ",")),
  sprintf("FastBAPS tip strip: %s cells; colored: %s", nrow(built$data[[fast_ring_layer_index]]), rendered_fast_strip),
  "Tree topology and supplied branch lengths were retained; only midpoint rooting changed display orientation.",
  "fastBAPS 7 and 8 are represented by multiple maximal local pure sectors; see Pseudomonas.fastBAPS_local_clade_sectors.csv."
)
writeLines(validation, file.path(out_dir, "validation.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat(paste(validation, collapse = "\n"), "\nOutput directory:", out_dir, "\n")
