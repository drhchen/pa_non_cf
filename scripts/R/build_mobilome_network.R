#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: build_mobilome_network.R DIST_TSV ELEMENTS_TSV METADATA_TSV OUTPUT_DIR",
    call. = FALSE
  )
}

required <- c("dplyr", "igraph", "readr", "tibble")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

dist_path <- args[[1]]
manifest_path <- args[[2]]
metadata_path <- args[[3]]
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dist_connection <- file(dist_path, open = "rb")
dist <- readr::read_tsv(
  dist_connection,
  col_names = c("query", "target", "mutation_distance", "p_value", "jaccard_raw"),
  col_types = "ccccc",
  progress = FALSE,
  trim_ws = TRUE
)
close(dist_connection)
read_tsv_portable <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection))
  readr::read_tsv(connection, show_col_types = FALSE, progress = FALSE)
}
manifest <- read_tsv_portable(manifest_path)
metadata <- read_tsv_portable(metadata_path)
if (!all(c("element_id", "sample", "element_type") %in% names(manifest))) {
  stop("Element manifest lacks element_id, sample, or element_type", call. = FALSE)
}
if (!all(c("sample", "phylogroup") %in% names(metadata))) {
  stop("Metadata lacks sample or phylogroup", call. = FALSE)
}
if (anyDuplicated(manifest$element_id) || anyDuplicated(metadata$sample)) {
  stop("Element IDs and metadata sample IDs must be unique", call. = FALSE)
}

normalize_element <- function(value) {
  value <- basename(value)
  sub("\\.(fasta|fna|fa)(\\.gz)?$", "", value, ignore.case = TRUE)
}
parse_jaccard <- function(value) {
  pieces <- strsplit(value, "/", fixed = TRUE)[[1]]
  if (length(pieces) != 2) {
    return(NA_real_)
  }
  numerator <- suppressWarnings(as.numeric(pieces[[1]]))
  denominator <- suppressWarnings(as.numeric(pieces[[2]]))
  if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  numerator / denominator
}

pairs <- dist |>
  dplyr::transmute(
    query = normalize_element(.data$query),
    target = normalize_element(.data$target),
    mutation_distance = suppressWarnings(as.numeric(.data$mutation_distance)),
    p_value = suppressWarnings(as.numeric(.data$p_value)),
    jaccard_index = vapply(.data$jaccard_raw, parse_jaccard, numeric(1))
  ) |>
  dplyr::filter(.data$query != .data$target, is.finite(.data$jaccard_index)) |>
  dplyr::mutate(
    from = pmin(.data$query, .data$target),
    to = pmax(.data$query, .data$target)
  ) |>
  dplyr::group_by(.data$from, .data$to) |>
  dplyr::summarise(
    jaccard_index = mean(.data$jaccard_index),
    mutation_distance = mean(.data$mutation_distance, na.rm = TRUE),
    p_value = max(.data$p_value, na.rm = TRUE),
    duplicate_rows_collapsed = dplyr::n(),
    .groups = "drop"
  )
if (nrow(pairs) == 0) {
  stop("No valid non-self BinDash comparisons were parsed", call. = FALSE)
}

unknown_elements <- setdiff(unique(c(pairs$from, pairs$to)), manifest$element_id)
if (length(unknown_elements) > 0) {
  stop(
    "BinDash names absent from the element manifest: ",
    paste(head(unknown_elements, 20), collapse = ", "),
    call. = FALSE
  )
}

mean_jaccard <- mean(pairs$jaccard_index)
edges <- pairs |>
  dplyr::filter(.data$jaccard_index >= mean_jaccard) |>
  dplyr::arrange(dplyr::desc(.data$jaccard_index), .data$from, .data$to)
if (nrow(edges) == 0) {
  stop("No edges remain at or above the mean Jaccard threshold", call. = FALSE)
}

network_elements <- unique(c(edges$from, edges$to))
vertices <- manifest |>
  dplyr::filter(.data$element_id %in% network_elements) |>
  dplyr::transmute(
    name = .data$element_id,
    sample = .data$sample,
    element_type = .data$element_type
  ) |>
  dplyr::left_join(metadata |> dplyr::select("sample", "phylogroup"), by = "sample") |>
  dplyr::mutate(phylogroup = dplyr::coalesce(as.character(.data$phylogroup), "Unknown"))

# Match the original graph definition: only elements in retained edges are nodes.
graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
igraph::E(graph)$weight <- edges$jaccard_index
components <- igraph::components(graph)
igraph::V(graph)$component <- components$membership
degree <- igraph::degree(graph)

density_value <- igraph::edge_density(graph, loops = FALSE)
clustering_value <- igraph::transitivity(graph, type = "average", isolates = "zero")
centralization_value <- igraph::centr_degree(
  graph,
  mode = "all",
  loops = FALSE,
  normalized = TRUE
)$centralization
heterogeneity_value <- if (mean(degree) == 0) NA_real_ else stats::sd(degree) / mean(degree)

# igraph <=1.x accepts the parameters reported in the manuscript. Newer
# versions removed them; retain a deterministic fallback and record sessionInfo.
set.seed(123)
layout <- tryCatch(
  igraph::layout_with_fr(
    graph,
    niter = 1500,
    area = 8000,
    repulserad = 1500,
    weights = igraph::E(graph)$weight
  ),
  error = function(error) {
    warning(
      "Installed igraph does not support area/repulserad; using the modern ",
      "Fruchterman-Reingold implementation. Match the publication igraph ",
      "version before regenerating the final figure."
    )
    igraph::layout_with_fr(graph, niter = 1500, weights = igraph::E(graph)$weight)
  }
)
layout <- layout * 1.3
layout_table <- tibble::tibble(
  element_id = igraph::V(graph)$name,
  x = layout[, 1],
  y = layout[, 2]
)

vertex_table <- igraph::as_data_frame(graph, what = "vertices") |>
  dplyr::mutate(degree = degree) |>
  dplyr::left_join(layout_table, by = c("name" = "element_id"))

readr::write_tsv(pairs, file.path(output_dir, "all_unique_pairs.tsv"))
readr::write_tsv(edges, file.path(output_dir, "network_edges.tsv"))
readr::write_tsv(vertex_table, file.path(output_dir, "network_vertices.tsv"))
readr::write_tsv(layout_table, file.path(output_dir, "network_layout.tsv"))
readr::write_tsv(
  tibble::tibble(
    metric = c(
      "mean_jaccard_threshold",
      "elements_total",
      "elements_in_network",
      "elements_without_retained_edge",
      "vertices",
      "edges",
      "isolated_vertices",
      "connected_components",
      "density",
      "average_clustering_coefficient",
      "degree_centralization",
      "degree_heterogeneity_cv"
    ),
    value = c(
      mean_jaccard,
      nrow(manifest),
      length(network_elements),
      nrow(manifest) - length(network_elements),
      igraph::vcount(graph),
      igraph::ecount(graph),
      sum(degree == 0),
      components$no,
      density_value,
      clustering_value,
      centralization_value,
      heterogeneity_value
    )
  ),
  file.path(output_dir, "network_metrics.tsv")
)

group_colors <- c(BAPS1 = "#FFD54F", BAPS2 = "#1E88E5", Unknown = "#BDBDBD")
vertex_colors <- unname(group_colors[vertex_table$phylogroup])
vertex_colors[is.na(vertex_colors)] <- group_colors[["Unknown"]]
vertex_shapes <- ifelse(
  vertex_table$element_type == "ICE",
  "square",
  ifelse(vertex_table$element_type == "IME", "circle", "diamond")
)
edge_widths <- 0.5 + 4 * (igraph::E(graph)$weight / max(igraph::E(graph)$weight))

grDevices::pdf(file.path(output_dir, "mobilome_network.pdf"), width = 10, height = 10)
plot(
  graph,
  layout = layout,
  vertex.size = 3,
  vertex.label = NA,
  vertex.color = vertex_colors,
  vertex.shape = vertex_shapes,
  vertex.frame.color = "white",
  edge.width = edge_widths,
  edge.color = grDevices::adjustcolor("#B0BEC5", alpha.f = 0.4),
  main = "ICE/IME similarity network"
)
legend(
  "topleft",
  legend = names(group_colors),
  col = unname(group_colors),
  pch = 16,
  pt.cex = 1.2,
  bty = "n",
  title = "Phylogroup"
)
legend(
  "topright",
  legend = c("ICE", "IME", "Unknown"),
  pch = c(15, 16, 18),
  col = "black",
  bty = "n",
  title = "Element type"
)
grDevices::dev.off()

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
