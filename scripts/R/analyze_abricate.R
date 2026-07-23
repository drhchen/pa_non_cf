#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3, 4)) {
  stop("Usage: analyze_abricate.R HITS_LONG METADATA OUTPUT_DIR [OVERRIDES]", call. = FALSE)
}

required <- c("dplyr", "ggplot2", "readr", "tibble", "tidyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

hits_path <- args[[1]]
metadata_path <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv_portable <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection))
  readr::read_tsv(connection, show_col_types = FALSE, progress = FALSE)
}
hits <- read_tsv_portable(hits_path)
metadata <- read_tsv_portable(metadata_path)
required_hit_columns <- c("sample", "database", "gene", "category")
required_metadata_columns <- c("sample", "phylogroup")
if (!all(required_hit_columns %in% names(hits))) {
  stop("Hit table is missing required columns", call. = FALSE)
}
if (!all(required_metadata_columns %in% names(metadata))) {
  stop("Metadata is missing required columns", call. = FALSE)
}
if (anyDuplicated(metadata$sample)) {
  stop("Metadata sample IDs are not unique", call. = FALSE)
}

metadata_two <- metadata |>
  dplyr::filter(.data$phylogroup %in% c("BAPS1", "BAPS2")) |>
  dplyr::mutate(phylogroup = factor(.data$phylogroup, levels = c("BAPS1", "BAPS2")))
if (!all(c("BAPS1", "BAPS2") %in% unique(as.character(metadata_two$phylogroup)))) {
  stop("Both BAPS1 and BAPS2 must be represented in metadata", call. = FALSE)
}
unknown_samples <- setdiff(unique(hits$sample), metadata$sample)
if (length(unknown_samples) > 0) {
  stop("Hits contain samples absent from metadata: ", paste(unknown_samples, collapse = ", "))
}

database_name <- if (nrow(hits) > 0) unique(hits$database) else "unknown"
if (length(database_name) != 1) {
  stop("The hit table contains more than one database", call. = FALSE)
}
if (length(args) == 4 && nzchar(args[[4]])) {
  overrides <- read_tsv_portable(args[[4]])
  required_override_columns <- c("database", "gene", "category", "reason", "source")
  if (!all(required_override_columns %in% names(overrides))) {
    stop("Override table is missing required columns", call. = FALSE)
  }
  relevant_overrides <- overrides |>
    dplyr::filter(tolower(.data$database) == tolower(database_name)) |>
    dplyr::select("gene", override_category = "category")
  if (anyDuplicated(relevant_overrides$gene)) {
    stop("Override table contains duplicate genes for this database", call. = FALSE)
  }
  hits <- hits |>
    dplyr::left_join(relevant_overrides, by = "gene") |>
    dplyr::mutate(category = dplyr::coalesce(.data$override_category, .data$category)) |>
    dplyr::select(-"override_category")
}

presence <- hits |>
  dplyr::distinct(.data$sample, .data$gene) |>
  dplyr::mutate(present = 1L)
genes <- sort(unique(hits$gene))
if (length(genes) > 0) {
  presence_matrix <- tidyr::expand_grid(sample = metadata$sample, gene = genes) |>
    dplyr::left_join(presence, by = c("sample", "gene")) |>
    dplyr::mutate(present = tidyr::replace_na(.data$present, 0L)) |>
    tidyr::pivot_wider(names_from = "gene", values_from = "present") |>
    dplyr::left_join(metadata, by = "sample")
} else {
  presence_matrix <- metadata
}
readr::write_tsv(presence_matrix, file.path(output_dir, "gene_presence_absence.tsv"))

test_gene <- function(gene_name) {
  sample_data <- metadata_two |>
    dplyr::select("sample", "phylogroup") |>
    dplyr::left_join(
      presence |> dplyr::filter(.data$gene == gene_name),
      by = "sample"
    ) |>
    dplyr::mutate(present = tidyr::replace_na(.data$present, 0L))
  contingency <- table(sample_data$phylogroup, factor(sample_data$present, levels = c(0, 1)))
  test <- stats::fisher.test(contingency)
  n_baps1 <- sum(sample_data$phylogroup == "BAPS1")
  n_baps2 <- sum(sample_data$phylogroup == "BAPS2")
  present_baps1 <- sum(sample_data$phylogroup == "BAPS1" & sample_data$present == 1)
  present_baps2 <- sum(sample_data$phylogroup == "BAPS2" & sample_data$present == 1)
  frequency_baps1 <- present_baps1 / n_baps1
  frequency_baps2 <- present_baps2 / n_baps2
  tibble::tibble(
    gene = gene_name,
    n_BAPS1 = n_baps1,
    n_BAPS2 = n_baps2,
    present_BAPS1 = present_baps1,
    present_BAPS2 = present_baps2,
    frequency_BAPS1 = frequency_baps1,
    frequency_BAPS2 = frequency_baps2,
    frequency_difference = frequency_baps2 - frequency_baps1,
    odds_ratio_BAPS2_vs_BAPS1 = unname(test$estimate),
    p_value = test$p.value
  )
}

if (length(genes) > 0) {
  gene_tests <- dplyr::bind_rows(lapply(genes, test_gene)) |>
    dplyr::mutate(q_value_bh = stats::p.adjust(.data$p_value, method = "BH")) |>
    dplyr::arrange(.data$q_value_bh, .data$p_value, .data$gene)
} else {
  gene_tests <- tibble::tibble()
}
readr::write_tsv(gene_tests, file.path(output_dir, "gene_fisher_tests.tsv"))

categories <- hits |>
  dplyr::filter(!is.na(.data$category), .data$category != "") |>
  tidyr::separate_longer_delim("category", delim = ";") |>
  dplyr::mutate(category = trimws(.data$category)) |>
  dplyr::filter(.data$category != "") |>
  dplyr::distinct(.data$sample, .data$gene, .data$category)

if (nrow(categories) > 0) {
  category_counts <- categories |>
    dplyr::count(.data$sample, .data$category, name = "n_unique_genes") |>
    tidyr::complete(
      sample = metadata$sample,
      category = sort(unique(categories$category)),
      fill = list(n_unique_genes = 0L)
    ) |>
    dplyr::left_join(metadata, by = "sample")
} else {
  category_counts <- tibble::tibble()
}
unique_counts <- presence |>
  dplyr::count(.data$sample, name = "n_unique_genes") |>
  tidyr::complete(sample = metadata$sample, fill = list(n_unique_genes = 0L)) |>
  dplyr::left_join(metadata, by = "sample") |>
  dplyr::mutate(category = "ALL_UNIQUE_GENES")
count_data <- dplyr::bind_rows(unique_counts, category_counts)
readr::write_tsv(count_data, file.path(output_dir, "category_counts.tsv"))

test_count <- function(category_name) {
  values <- count_data |>
    dplyr::filter(
      .data$category == category_name,
      .data$phylogroup %in% c("BAPS1", "BAPS2")
    ) |>
    dplyr::mutate(phylogroup = factor(.data$phylogroup, levels = c("BAPS1", "BAPS2")))
  test <- stats::wilcox.test(n_unique_genes ~ phylogroup, data = values, exact = FALSE)
  median_baps1 <- stats::median(values$n_unique_genes[values$phylogroup == "BAPS1"])
  median_baps2 <- stats::median(values$n_unique_genes[values$phylogroup == "BAPS2"])
  tibble::tibble(
    category = category_name,
    median_BAPS1 = median_baps1,
    median_BAPS2 = median_baps2,
    median_difference_BAPS2_minus_BAPS1 = median_baps2 - median_baps1,
    p_value = test$p.value
  )
}
count_tests <- dplyr::bind_rows(lapply(sort(unique(count_data$category)), test_count)) |>
  dplyr::mutate(q_value_bh = stats::p.adjust(.data$p_value, method = "BH")) |>
  dplyr::arrange(.data$q_value_bh, .data$category)
readr::write_tsv(count_tests, file.path(output_dir, "category_wilcoxon_tests.tsv"))

plot_data <- count_data |>
  dplyr::filter(.data$phylogroup %in% c("BAPS1", "BAPS2")) |>
  dplyr::mutate(phylogroup = factor(.data$phylogroup, levels = c("BAPS1", "BAPS2")))
plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = .data$phylogroup, y = .data$n_unique_genes, fill = .data$phylogroup)
) +
  ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.75) +
  ggplot2::geom_jitter(width = 0.12, height = 0, size = 0.7, alpha = 0.5) +
  ggplot2::facet_wrap(~category, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(BAPS1 = "#E69F00", BAPS2 = "#56B4E9")) +
  ggplot2::labs(x = NULL, y = "Unique detected genes") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave(file.path(output_dir, "category_counts.pdf"), plot, width = 10, height = 7)

writeLines(
  c(
    paste0("database\t", database_name),
    paste0("samples_in_metadata\t", nrow(metadata)),
    paste0("samples_in_two_lineage_test\t", nrow(metadata_two)),
    paste0("raw_hits\t", nrow(hits)),
    paste0("unique_sample_gene_pairs\t", nrow(presence)),
    paste0("unique_genes\t", length(genes)),
    "count_unit\tunique genes per sample; repeated loci/hits are collapsed"
  ),
  file.path(output_dir, "analysis_summary.tsv")
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
