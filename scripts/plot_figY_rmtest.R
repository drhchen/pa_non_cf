# Fig. Y — lineage-resolved recombination landscape + r/m permutation test
# Inputs  : gubbins/hotspot_reproducibility/figY_*.csv, rm_test_input.csv
# Outputs : gubbins/figures/FigY_recombination_landscape.{png,pdf}
#           gubbins/hotspot_reproducibility/rm_permutation_test.txt

Sys.setenv(LANG = "en")
library(tidyverse)
source("results/baps_palette.R")

root <- normalizePath(getwd(), mustWork=TRUE)
setwd(root)

out_dir <- "results/3.2_recombination_gubbins/figures"
dir.create(out_dir, showWarnings = FALSE)
res_dir <- "/tmp/figy_inputs2"

## ---- load data ---------------------------------------------------------------
cov <- read_csv(file.path(res_dir, "figY_coverage_binned.csv"),
                show_col_types = FALSE) %>%
  transmute(bin_start, bin_end, mid = (bin_start + bin_end) / 2,
            baps7 = coverage_baps7, baps8 = coverage_baps8) %>%
  pivot_longer(c(baps7, baps8), names_to = "lineage", values_to = "coverage")

hot <- read_csv(file.path(res_dir, "figY_hotspot_tracks.csv"),
                show_col_types = FALSE) %>%
  mutate(lineage = factor(lineage, levels = c("baps7", "baps8")))

recurrent_spans <- hot %>%
  group_by(lineage, class, start, end) %>% slice(1) %>% ungroup()

# labels for the five lineage-independent hotspots
rec_labels <- recurrent_spans %>%
  filter(class == "lineage-independent") %>%
  group_by(start, end) %>% slice(1) %>% ungroup() %>%
  mutate(mid = (start + end) / 2,
         label = c("cat / Pf1 phage",
                   "pyoS5 / cupC",
                   "fabF2 / pmpM",
                   "cupB / bfmRS / ampR-ampC",
                   "nadC / pilA-pilB"))

lin_colors <- c(baps7 = BAPS_COLORS[["BAPS7"]], baps8 = BAPS_COLORS[["BAPS8"]])

# gold shading must appear in BOTH panels: duplicate per facet level
gold_bands <- bind_rows(
  mutate(rec_labels, lineage = factor("baps7", levels = c("baps7", "baps8"))),
  mutate(rec_labels, lineage = factor("baps8", levels = c("baps7", "baps8"))))

p <- ggplot() +
  # lineage-specific detections: light marks in their own track
  geom_rect(data = filter(recurrent_spans, class == "lineage-specific"),
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
                fill = "lineage-specific detection"), alpha = 0.18) +
  # lineage-independent hotspots: strong shading spanning both tracks
  geom_rect(data = gold_bands,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "goldenrod1", alpha = 0.25, inherit.aes = FALSE) +
  # coverage tracks
  geom_col(data = cov, aes(x = mid, y = coverage, group = lineage,
                           fill = lineage),
           width = 10000, colour = "white", linewidth = 0.05) +
  facet_wrap(~lineage, ncol = 1, strip.position = "right") +
  scale_fill_manual(values = c(lin_colors,
                               "lineage-specific detection" = "grey55"),
                    breaks = c("baps7", "baps8", "lineage-specific detection"),
                    name = "Lineage",
                    labels = c("baps7", "baps8", "Lineage-specific hotspot")) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6, "M"),
                     expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.32))) +
  coord_cartesian(clip = "off") +
  labs(x = "PAO1 genomic position (bp)",
       y = "Recombination events per 10 kb bin",
       title = "Lineage-resolved recombination landscape",
       subtitle = "Gold shading: lineage-independent hotspots (detected in both lineages); grey shading: lineage-specific detections") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold"),
        legend.position = "top",
        plot.title.position = "plot")

# hotspot labels anchored at the left edge of each band, angled, inside panel
p <- p + geom_text(data = mutate(rec_labels, lineage = "baps8"),
                   aes(x = start, y = Inf, label = label),
                   inherit.aes = FALSE, vjust = 1.1,
                   size = 2.6, angle = 32, hjust = 0, colour = "grey15",
                   check_overlap = FALSE)

ggsave(file.path(out_dir, "FigY_recombination_landscape.png"),
       p, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "FigY_recombination_landscape.pdf"),
       p, width = 10, height = 6.5)

## ---- r/m permutation test -----------------------------------------------------
rm <- read_csv(file.path(res_dir, "rm_test_input.csv"), show_col_types = FALSE)
v7 <- rm %>% filter(lineage == "baps7") %>% pull(r_over_m)
v8 <- rm %>% filter(lineage == "baps8") %>% pull(r_over_m)

observed_median <- median(v8) - median(v7)
observed_mean <- mean(v8) - mean(v7)
set.seed(20261109)
B <- 10000
perm_median <- numeric(B)
perm_mean <- numeric(B)
pooled <- c(v7, v8); n7 <- length(v7)
for (b in seq_len(B)) {
  idx <- sample(length(pooled))
  perm_median[b] <- median(pooled[idx[(n7 + 1):length(pooled)]]) -
    median(pooled[idx[1:n7]])
  perm_mean[b] <- mean(pooled[idx[(n7 + 1):length(pooled)]]) -
    mean(pooled[idx[1:n7]])
}
p_perm_median <- mean(abs(perm_median) >= abs(observed_median))
p_perm_mean <- mean(abs(perm_mean) >= abs(observed_mean))

# Hodges-Lehmann shift estimate + Mann-Whitney
hl <- median(as.vector(outer(v8, v7, "-")))
mw <- suppressWarnings(wilcox.test(v8, v7))

# effect size: Cliff's delta
cliff <- (sum(outer(v8, v7, ">")) - sum(outer(v8, v7, "<"))) / (length(v8) * length(v7))

report <- sprintf(
"Per-branch r/m comparison (tip branches only)
================================================
baps7: n = %d tips, median = %.4f, IQR = %.4f-%.4f, mean = %.4f
baps8: n = %d tips, median = %.4f, IQR = %.4f-%.4f, mean = %.4f

Median difference (baps8 - baps7): %.4f
Mean difference (baps8 - baps7):   %.4f
Permutation test, median-diff statistic (two-sided, B = %d): p = %.4f
Permutation test, mean-diff statistic   (two-sided, B = %d): p = %.4f
Hodges-Lehmann shift: %.4f
Mann-Whitney U: p = %.3g
Cliff's delta: %.3f (%.0f%% of baps8 tips above the baps7 median)

Caveat: tip branches within a lineage share ancestry, so branch-level
values are not fully independent; the permutation p-values are
approximate and the effect sizes should be reported alongside them.
",
length(v7), median(v7), quantile(v7, .25), quantile(v7, .75), mean(v7),
length(v8), median(v8), quantile(v8, .25), quantile(v8, .75), mean(v8),
observed_median, observed_mean, B, p_perm_median, B, p_perm_mean,
hl, mw$p.value, cliff,
100 * mean(v8 > median(v7)))

write_lines(report, file.path(res_dir, "rm_permutation_test.txt"))
cat(report)
