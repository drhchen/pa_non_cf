#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: infer_population_structure.R ALIGNMENT OUTPUT_DIR MAX_DEPTH MAX_POPS",
    call. = FALSE
  )
}

required <- c("readr", "rhierbaps", "tibble")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

alignment_path <- args[[1]]
output_dir <- args[[2]]
max_depth <- as.integer(args[[3]])
max_pops <- as.integer(args[[4]])
if (!file.exists(alignment_path) || is.na(max_depth) || is.na(max_pops)) {
  stop("Invalid input path or numeric argument", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

snp_matrix <- rhierbaps::load_fasta(alignment_path)
if (nrow(snp_matrix) < 2 || ncol(snp_matrix) < 1) {
  stop("Alignment must contain at least two samples and one site", call. = FALSE)
}

set.seed(123)
fit <- rhierbaps::hierBAPS(
  snp_matrix,
  max.depth = max_depth,
  n.pops = max_pops,
  quiet = TRUE
)
partition <- as.data.frame(fit$partition.df, check.names = FALSE)
if (nrow(partition) != nrow(snp_matrix)) {
  stop("rhierBAPS returned a partition with an unexpected number of samples", call. = FALSE)
}
if (!any(names(partition) %in% c("sample", "Isolate", "isolate"))) {
  partition <- tibble::rownames_to_column(partition, var = "sample")
}

readr::write_tsv(partition, file.path(output_dir, "hierbaps_partitions.tsv"))
saveRDS(fit, file.path(output_dir, "hierbaps_fit.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
