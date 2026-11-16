#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
source("results/baps_palette.R")
})

input_dir <- "results/3.5_mobilome/network_analysis_all_elements_jaccard_gt_0.8"
# This run uses the pre-filtered all-elements files supplied by the user.
out_dir <- input_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
threshold <- 0.8

# FastBAPS file is comma-separated in this project; support either delimiter.
fb_file <- "baps/fastbaps_result.csv"
fb <- read.csv(fb_file, stringsAsFactors = FALSE, check.names = FALSE)
names(fb)[1:2] <- c("sample", "group")
group_key <- setNames(as.character(fb$group), as.character(fb$sample))

extract_sample <- function(x) {
  x <- sub("^\\./", "", x)
  y <- x
  out <- rep(NA_character_, length(y))
  m <- regexec("^((ERR|SRR)[0-9]+|LHL[-_][0-9]+|(GCA|GCF)_[0-9.]+|PAO1)(_|-|$)", y, perl = TRUE)
  z <- regmatches(y, m)
  for (i in seq_along(z)) if (length(z[[i]]) >= 2) out[i] <- z[[i]][2]
  out
}

files <- sort(list.files(input_dir, pattern = "_dist_all_elements_jaccard_gt_0\\.8\\.txt$", full.names = TRUE))
if (!length(files)) stop("No *_dist.txt files found in ", input_dir)

all_summary <- list()
for (f in files) {
  feature <- sub("_dist_all_elements_jaccard_gt_0\\.8\\.txt$", "", basename(f))
  d <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("query", "target", "mutation_distance", "p_value", "jaccard_raw"),
                  quote = "", comment.char = "")
  parts <- strsplit(d$jaccard_raw, "/", fixed = TRUE)
  d$jaccard_index <- vapply(parts, function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
  edges <- d[d$jaccard_index > threshold, c("query", "target", "mutation_distance", "p_value", "jaccard_index")]
  # Remove self-loops and guard against accidental duplicate pairs.
  edges <- edges[edges$query != edges$target, ]
  if (nrow(edges)) {
    key <- paste(pmin(edges$query, edges$target), pmax(edges$query, edges$target), sep = "\r")
    edges <- edges[!duplicated(key), ]
  }
  nodes <- sort(unique(c(d$query, d$target)))
  g <- graph_from_data_frame(edges[, c("query", "target", "jaccard_index", "mutation_distance"), drop = FALSE],
                             directed = FALSE, vertices = data.frame(name = nodes))
  if (ecount(g)) {
    E(g)$weight <- edges$jaccard_index
    E(g)$mutation_distance <- edges$mutation_distance
  }
  deg <- degree(g)
  density_val <- edge_density(g, loops = FALSE)
  clustering_coeff <- if (vcount(g) > 2 && ecount(g) > 0) transitivity(g, type = "average", isolates = "zero") else 0
  centralization_val <- if (vcount(g) > 1) centr_degree(g, normalized = TRUE)$centralization else 0
  heterogeneity_val <- if (length(deg) && mean(deg) > 0) sd(deg) / mean(deg) else 0
  cl <- components(g)
  V(g)$sample <- extract_sample(V(g)$name)
  V(g)$group <- unname(group_key[V(g)$sample]); V(g)$group[is.na(V(g)$group)] <- "Unknown"
  V(g)$cluster_id <- cl$membership
  component_sizes <- sort(as.integer(cl$csize), decreasing = TRUE)
  comp_df <- data.frame(feature = feature, cluster_id = seq_along(cl$csize), size = as.integer(cl$csize))
  write.table(comp_df, file.path(out_dir, paste0(feature, "_components.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(edges, file.path(out_dir, paste0(feature, "_edges_jaccard_gt_0.8.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(node = V(g)$name, sample = V(g)$sample, group = V(g)$group, cluster_id = V(g)$cluster_id, degree = deg),
              file.path(out_dir, paste0(feature, "_nodes.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  # Cluster-composition visualization: each bar is one connected component,
  # and colors show the number of elements from each FastBAPS group. This is
  # more readable than drawing all edges for dense networks.
  nd <- data.frame(cluster_id = V(g)$cluster_id, group = V(g)$group)
  tab <- as.data.frame(table(nd$cluster_id, nd$group), stringsAsFactors = FALSE)
  names(tab) <- c("cluster_id", "group", "count")
  tab$cluster_id <- as.integer(as.character(tab$cluster_id))
  cs <- sort(table(nd$cluster_id), decreasing = TRUE)
  top_ids <- as.integer(names(cs)[seq_len(min(30, length(cs)))])
  tab <- tab[tab$cluster_id %in% top_ids & tab$count > 0, ]
  tab$cluster_id <- factor(tab$cluster_id, levels = rev(top_ids))
  groups <- sort(unique(nd$group)); group_names <- ifelse(grepl("^BAPS", groups), groups, paste0("BAPS", groups)); cols <- setNames(BAPS_COLORS[group_names], groups); cols[is.na(cols)] <- "#BDBDBD"
  png(file.path(out_dir, paste0(feature, "_cluster_composition.png")), width = 3600, height = 2200, res = 260)
  par(mar = c(8, 5, 4, 2))
  mat <- xtabs(count ~ cluster_id + group, tab)
  mat <- mat[rev(rownames(mat)), , drop = FALSE]
  barplot(t(mat), beside = FALSE, col = cols[colnames(mat)], las = 2,
          ylab = "Number of elements", xlab = "Connected component (ranked by size)",
          main = paste0(feature, ": FastBAPS composition of largest connected components"),
          legend.text = colnames(mat), args.legend = list(x = "topright", bty = "n", cex = 0.8))
  dev.off()
  pdf(file.path(out_dir, paste0(feature, "_cluster_composition.pdf")), width = 14, height = 9, useDingbats = FALSE)
  par(mar = c(8, 5, 4, 2))
  barplot(t(mat), beside = FALSE, col = cols[colnames(mat)], las = 2,
          ylab = "Number of elements", xlab = "Connected component (ranked by size)",
          main = paste0(feature, ": FastBAPS composition of largest connected components"),
          legend.text = colnames(mat), args.legend = list(x = "topright", bty = "n", cex = 0.8))
  dev.off()

  # Reproducible force-directed plot. Components are shaded and FastBAPS groups colored.
  set.seed(123)
  # Expand the coordinates substantially so large, dense components remain legible.
  # The larger canvas below preserves this spacing in the exported PNG.
  lay <- if (vcount(g) > 1) layout_with_fr(g, niter = 3000) else matrix(c(0, 0), ncol = 2)
  lay <- lay * if (vcount(g) > 300) 3.5 else if (vcount(g) > 100) 2.5 else 1.8
  # Force very large components onto expanded rings. This prevents the dense core
  # from collapsing into an unreadable ball of overlapping vertices.
  memb <- components(g)$membership
  for (cc in unique(memb)) {
    ii <- which(memb == cc); nn <- length(ii)
    if (nn > 50) {
      ctr <- colMeans(lay[ii, , drop = FALSE])
      ring <- layout_in_circle(induced_subgraph(g, vids = ii))
      ring <- ring * (2.5 + sqrt(nn) / 2)
      lay[ii, ] <- sweep(ring, 2, ctr, "+")
    }
  }
  grp_levels <- sort(unique(V(g)$group)); pal <- setNames(grDevices::hcl.colors(max(1, length(grp_levels)), "Set 2"), grp_levels)
  vcol <- unname(pal[V(g)$group]); vcol[V(g)$group == "Unknown"] <- "grey75"
  png(file.path(out_dir, paste0(feature, "_network.png")), width = 4000, height = 3200, res = 260)
  plot(g, layout = lay, vertex.size = if (vcount(g) > 300) 2 else 2.5, vertex.label = NA, vertex.color = vcol,
       vertex.frame.color = "white", vertex.frame.width = 0.25,
       edge.width = if (ecount(g)) 0.4 + 2 * E(g)$weight else 0.5,
       edge.color = adjustcolor("#607D8B", alpha.f = 0.18), main = paste0(feature, " network (Jaccard > ", threshold, ")"))
  legend("topleft", legend = names(pal), col = pal, pch = 16, pt.cex = 1.1, bty = "n", title = "FastBAPS cluster")
  dev.off()
  pdf(file.path(out_dir, paste0(feature, "_network.pdf")), width = 14, height = 11, useDingbats = FALSE)
  plot(g, layout = lay, vertex.size = if (vcount(g) > 300) 2 else 2.5, vertex.label = NA, vertex.color = vcol,
       vertex.frame.color = "white", vertex.frame.width = 0.25,
       edge.width = if (ecount(g)) 0.4 + 2 * E(g)$weight else 0.5,
       edge.color = adjustcolor("#607D8B", alpha.f = 0.18), main = paste0(feature, " network (Jaccard > ", threshold, ")"))
  legend("topleft", legend = names(pal), col = pal, pch = 16, pt.cex = 1.1, bty = "n", title = "FastBAPS cluster")
  dev.off()

  # Dedicated high-resolution view of the largest connected component.
  # This is essential for dense networks where the full graph compresses the core.
  largest_id <- which.max(cl$csize)
  core <- induced_subgraph(g, vids = which(cl$membership == largest_id))
  if (vcount(core) > 1) {
    set.seed(456)
    core_lay <- layout_with_fr(core, niter = 5000)
    core_lay <- core_lay * 8
    core_groups <- V(core)$group
    core_col <- unname(pal[core_groups]); core_col[is.na(core_col)] <- "grey75"
    png(file.path(out_dir, paste0(feature, "_largest_component_zoom.png")), width = 5000, height = 5000, res = 300)
    plot(core, layout = core_lay, vertex.size = if (vcount(core) > 300) 3 else 5,
         vertex.label = if (vcount(core) <= 80) V(core)$sample else NA,
         vertex.label.cex = 0.45, vertex.label.color = "grey20",
         vertex.color = core_col, vertex.frame.color = "white", vertex.frame.width = 0.3,
         edge.width = if (ecount(core)) 0.5 + 2.5 * E(core)$weight else 0.5,
         edge.color = adjustcolor("#455A64", alpha.f = 0.28),
         main = paste0(feature, " largest component (", vcount(core), " nodes; Jaccard > ", threshold, ")"))
    legend("topleft", legend = names(pal), col = pal, pch = 16, pt.cex = 1.2, bty = "n", title = "FastBAPS cluster")
    dev.off()
    pdf(file.path(out_dir, paste0(feature, "_largest_component_zoom.pdf")), width = 12, height = 12, useDingbats = FALSE)
    plot(core, layout = core_lay, vertex.size = if (vcount(core) > 300) 3 else 5,
         vertex.label = if (vcount(core) <= 80) V(core)$sample else NA,
         vertex.label.cex = 0.45, vertex.label.color = "grey20",
         vertex.color = core_col, vertex.frame.color = "white", vertex.frame.width = 0.3,
         edge.width = if (ecount(core)) 0.5 + 2.5 * E(core)$weight else 0.5,
         edge.color = adjustcolor("#455A64", alpha.f = 0.28),
         main = paste0(feature, " largest component (", vcount(core), " nodes; Jaccard > ", threshold, ")"))
    legend("topleft", legend = names(pal), col = pal, pch = 16, pt.cex = 1.2, bty = "n", title = "FastBAPS cluster")
    dev.off()
    # Optional sparse companion view is disabled; the requested output retains
    # the complete network edge set.
    if (FALSE && ecount(core) > 0) {
      keep <- logical(ecount(core)); inc <- incident_edges(core, V(core), mode = "all")
      for (ii in seq_len(vcount(core))) {
        ee <- inc[[ii]]; if (length(ee)) keep[ee[order(E(core)[ee]$weight, decreasing = TRUE)[seq_len(min(3, length(ee)))]]] <- TRUE
      }
      sparse <- subgraph.edges(core, which(keep), delete.vertices = FALSE)
      slay <- layout_with_fr(sparse, niter = 5000) * 7
      png(file.path(out_dir, paste0(feature, "_largest_component_top3_edges.png")), width = 5000, height = 5000, res = 300)
      scol <- unname(pal[V(sparse)$group]); scol[is.na(scol)] <- "grey75"
      plot(sparse, layout = slay, vertex.size = 5, vertex.label = NA, vertex.color = scol,
           vertex.frame.color = "white", edge.width = 1.2, edge.color = adjustcolor("#455A64", alpha.f = 0.45),
           main = paste0(feature, " largest component (top 3 edges per node)"))
      legend("topleft", legend = names(pal), col = pal, pch = 16, pt.cex = 1.2, bty = "n", title = "FastBAPS cluster")
      dev.off()
      pdf(file.path(out_dir, paste0(feature, "_largest_component_top3_edges.pdf")), width = 12, height = 12, useDingbats = FALSE)
      plot(sparse, layout = slay, vertex.size = 5, vertex.label = NA, vertex.color = scol,
           vertex.frame.color = "white", edge.width = 1.2, edge.color = adjustcolor("#455A64", alpha.f = 0.45),
           main = paste0(feature, " largest component (top 3 edges per node)"))
      dev.off()
    }
  }
  all_summary[[feature]] <- data.frame(feature = feature, total_pairs = nrow(d), nodes = vcount(g), edges_jaccard_gt_0.8 = ecount(g),
    components = components(g)$no, largest_component = ifelse(length(component_sizes), component_sizes[1], 0),
    network_density = density_val, clustering_coefficient = clustering_coeff,
    centralization = centralization_val, heterogeneity = heterogeneity_val)
}
summary_df <- do.call(rbind, all_summary)
write.table(summary_df, file.path(out_dir, "network_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(summary_df, row.names = FALSE)
