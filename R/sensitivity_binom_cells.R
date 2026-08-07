# sensitivity_binom_cells.R
#
# Adds binomial cell sampling to the main single-chromosome simulation.
# This addresses the paper's Limitations note that incorporating Binomial
# sampling of 5-10 biopsied cells raises the false-normal rate from
# 71.5% to 84.2% at M = 15%, 10× coverage.
#
# The original simulate_cnv() treats the biopsied mosaicism fraction as
# exactly equal to the embryo-wide fraction M. In reality, n_cells cells
# are sampled from the trophectoderm, and the number carrying the CNV is
# Binomial(n_cells, M), making the measured fraction random. This adds a
# layer of sampling noise *on top of* sequencing noise.
#
# This script produces:
#   results/binom_sensitivity.csv    — full results grid
#   figures/supp_binom_cells.png     — comparison figure (if ggplot2 present)
#
# The main conclusions of the paper are not weakened by this addition:
# false-normal rates increase, which makes the argument for probabilistic
# reporting *stronger*, not weaker.
#
# Author: Anere Cye

source("R/simulate_cnv_binom.R")
source("R/classify_embryo.R")

set.seed(2025)

# ── Parameters ─────────────────────────────────────────────────────────────
n_cells_values <- c(NA, 5, 8, 10)   # NA = original model (no binom sampling)
m_levels       <- c(0, 0.10, 0.15, 0.20, 0.25, 0.30)
coverages      <- c(5, 10, 20, 50)
low_thr        <- 0.15
high_thr       <- 0.30
n_replicates   <- 5000

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ── Wilson CI helper ────────────────────────────────────────────────────────
wilson_ci <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  ci <- prop.test(k, n, conf.level = conf, correct = FALSE)$conf.int
  as.numeric(ci)
}

# ── Run simulation ──────────────────────────────────────────────────────────
rows <- list()

for (nc in n_cells_values) {
  nc_label <- if (is.na(nc)) "fixed (original)" else as.character(nc)
  cat("n_cells =", nc_label, "...\n")

  for (cov in coverages) {
    for (m in m_levels) {
      log2r <- replicate(n_replicates,
                         simulate_cnv_binom(
                           true_mosaicism = m,
                           coverage       = cov,
                           n_cells        = nc
                         ))
      cls    <- classify_embryo(log2r, low_thr, high_thr)
      n_norm <- sum(cls == "normal")
      ci     <- wilson_ci(n_norm, n_replicates)

      rows[[length(rows) + 1]] <- data.frame(
        n_cells        = nc_label,
        coverage       = cov,
        true_mosaicism = m,
        prop_normal    = n_norm / n_replicates,
        ci_lower       = ci[1],
        ci_upper       = ci[2]
      )
    }
  }
}

results_binom <- do.call(rbind, rows)
write.csv(results_binom, "results/binom_sensitivity.csv", row.names = FALSE)
cat("Saved: results/binom_sensitivity.csv\n")

# ── Optional figure ─────────────────────────────────────────────────────────
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  plot_data <- results_binom[
    results_binom$coverage == 10 &
    results_binom$true_mosaicism %in% c(0, 0.10, 0.15, 0.20, 0.25, 0.30),
  ]
  plot_data$m_label <- paste0(plot_data$true_mosaicism * 100, "% mosaicism")

  # Ordered factor so legend is logical
  nc_order <- c("fixed (original)", "5", "8", "10")
  plot_data$n_cells <- factor(plot_data$n_cells,
                              levels = nc_order[nc_order %in% plot_data$n_cells])

  p <- ggplot(plot_data,
              aes(x = true_mosaicism * 100, y = prop_normal,
                  colour = n_cells, group = n_cells)) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.8, alpha = 0.5) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    scale_colour_brewer(palette = "Set1",
                        labels = c("fixed (original)" = "No biopsy sampling (original)",
                                   "5" = "5 cells", "8" = "8 cells", "10" = "10 cells")) +
    labs(
      x       = "True mosaicism level (%)",
      y       = "Proportion classified as normal",
      colour  = "Biopsy size",
      title   = "Supplementary: effect of binomial biopsy cell sampling",
      subtitle = "Coverage = 10×; thresholds 0.15 / 0.30; n = 5,000 per cell"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")

  ggsave("figures/supp_binom_cells.png", p,
         width = 8, height = 5, dpi = 300)
  cat("Saved: figures/supp_binom_cells.png\n")
} else {
  cat("ggplot2 not available — skipping figure.\n")
}

# ── Console summary (replicate Table 1 comparison) ─────────────────────────
cat("\n=== False-normal rate at M = 15%, coverage = 10× ===\n")
sub <- results_binom[
  results_binom$true_mosaicism == 0.15 &
  results_binom$coverage == 10,
  c("n_cells", "prop_normal", "ci_lower", "ci_upper")
]
sub$prop_normal_pct <- paste0(round(sub$prop_normal * 100, 1), "%")
sub$ci              <- paste0("[",
                               round(sub$ci_lower * 100, 1), "–",
                               round(sub$ci_upper * 100, 1), "%]")
print(sub[, c("n_cells", "prop_normal_pct", "ci")], row.names = FALSE)
