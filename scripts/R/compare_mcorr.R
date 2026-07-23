#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: compare_mcorr.R BAPS1_FIT BAPS2_FIT OUTPUT_PREFIX", call. = FALSE)
}

required <- c("dplyr", "ggplot2", "readr", "tibble", "tidyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

read_fit <- function(path, lineage) {
  connection <- file(path, open = "rb")
  on.exit(close(connection))
  data <- readr::read_csv(connection, show_col_types = FALSE, progress = FALSE)
  required_columns <- c("group", "ratio", "fbar", "c", "d_sample")
  if (!all(required_columns %in% names(data))) {
    stop(path, " is missing required mcorr columns", call. = FALSE)
  }
  data |>
    dplyr::filter(as.character(.data$group) != "all") |>
    dplyr::transmute(
      phylogroup = lineage,
      ratio = as.numeric(.data$ratio),
      fbar = as.numeric(.data$fbar),
      coverage = as.numeric(.data$c),
      diversity = as.numeric(.data$d_sample)
    )
}

data <- dplyr::bind_rows(read_fit(args[[1]], "BAPS1"), read_fit(args[[2]], "BAPS2")) |>
  dplyr::mutate(phylogroup = factor(.data$phylogroup, levels = c("BAPS1", "BAPS2")))
long <- data |>
  tidyr::pivot_longer(
    cols = c("diversity", "ratio", "fbar", "coverage"),
    names_to = "parameter",
    values_to = "value"
  ) |>
  dplyr::mutate(
    parameter = factor(
      .data$parameter,
      levels = c("diversity", "ratio", "fbar", "coverage"),
      labels = c(
        "d (sample diversity)",
        "gamma/mu (mcorr ratio)",
        "Mean recombination fragment size (bp)",
        "Recombination coverage"
      )
    )
  )

test_parameter <- function(parameter_name) {
  subset <- long |> dplyr::filter(.data$parameter == parameter_name, is.finite(.data$value))
  test <- stats::wilcox.test(value ~ phylogroup, data = subset, exact = FALSE)
  median_baps1 <- stats::median(subset$value[subset$phylogroup == "BAPS1"])
  median_baps2 <- stats::median(subset$value[subset$phylogroup == "BAPS2"])
  tibble::tibble(
    parameter = as.character(parameter_name),
    n_BAPS1 = sum(subset$phylogroup == "BAPS1"),
    n_BAPS2 = sum(subset$phylogroup == "BAPS2"),
    median_BAPS1 = median_baps1,
    median_BAPS2 = median_baps2,
    median_difference_BAPS2_minus_BAPS1 = median_baps2 - median_baps1,
    p_value = test$p.value
  )
}
tests <- dplyr::bind_rows(lapply(levels(long$parameter), test_parameter)) |>
  dplyr::mutate(q_value_bh = stats::p.adjust(.data$p_value, method = "BH"))

output_prefix <- args[[3]]
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(long, paste0(output_prefix, "_long.tsv"))
readr::write_tsv(tests, paste0(output_prefix, "_wilcoxon.tsv"))

labels <- tests |>
  dplyr::mutate(label = paste0("BH q = ", format(.data$q_value_bh, digits = 2, scientific = TRUE)))
plot <- ggplot2::ggplot(
  long,
  ggplot2::aes(x = .data$phylogroup, y = .data$value, fill = .data$phylogroup)
) +
  ggplot2::geom_violin(alpha = 0.55, scale = "width", trim = FALSE) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
  ggplot2::geom_jitter(width = 0.12, size = 0.65, alpha = 0.5) +
  ggplot2::geom_text(
    data = labels,
    ggplot2::aes(x = 1.5, y = Inf, label = .data$label),
    inherit.aes = FALSE,
    vjust = 1.3,
    size = 3
  ) +
  ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 2) +
  ggplot2::scale_fill_manual(values = c(BAPS1 = "#4E79A7", BAPS2 = "#F28E2B")) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave(paste0(output_prefix, ".pdf"), plot, width = 9, height = 7)
writeLines(capture.output(sessionInfo()), paste0(output_prefix, "_sessionInfo.txt"))
