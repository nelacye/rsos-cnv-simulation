# sensitivity_eta.R
#
# Sensitivity analysis for the PCR noise multiplier eta (pcr_duplicates).
#
# In the main simulation, eta is fixed at 1.2 (Sabina & Leamon, 2015).
# This parameter amplifies baseline noise and captures whole-genome
# amplification variance, but it varies across protocols and platforms.
# sigma_0 was already varied in Supplementary Figure S1; eta was not.
# This script fills that gap.
#
# Design:
#   - eta in {1.0, 1.2, 1.5, 2.0}
#   - M in {0, 0.10, 0.15, 0.20} (key mosaicism levels)
#   - Coverage in {5, 10, 20, 50}
#   - Thresholds: 0.15 / 0.30 (paper defaults)
#   - n_replicates = 5000 per cell (matching main simulation)
#
# Output: results/sensitivity_eta.csv
#         figures/supp_eta_sensitivity.png
#
# Author: Anere Cye

source("R/simulate_cnv_binom.R")   # uses simulate_cnv_binom()
source("R/classify_embryo.R")

set.seed(2025)

# ── Parameters ─────────────────────────────────────────────────────────────
eta_values    <- c(1.0, 1.2, 1.5, 2.0)
m_levels      <- c(0, 0.10, 0.15, 0.20)
coverages     <- c(5, 10, 20, 50)
low_thr       <- 0.15
high_thr      <- 0.30
n_replicates  <- 5000
n_cells_fixed <- NA    # NA = no binomial sampling; mirrors original model

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

for (eta in eta_values) {
  for (cov in coverages) {
    for (m in m_levels) {
      log2r  <- replicate(n_replicates,
                          simulate_cnv_binom(
                            true_mosaicism = m,
                            coverage       = cov,
                            pcr_duplicates = eta,
                            n_cells        = n_cells_fixed
                          ))
      cls    <- classify_embryo(log2r, low_thr, high_thr)
      n_norm <- sum(cls == "normal")
      ci     <- wilson_ci(n_norm, n_replicates)

      rows[[length(rows) + 1]] <- data.frame(
        eta            = eta,
        coverage       = cov,
        true_mosaicism = m,
        prop_normal    = n_norm / n_replicates,
        ci_lower       = ci[1],
        ci_upper       = ci[2]
      )
    }
  }
}

results_eta <- do.call(rbind, rows)
write.csv(results_eta, "results/sensitivity_eta.csv", row.names = FALSE)
cat("Saved: results/sensitivity_eta.csv\n")

# ── Optional figure (requires ggplot2) ─────────────────────────────────────
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  results_eta$eta_label <- paste0("η = ", results_eta$eta)
  results_eta$m_label   <- paste0("M = ", results_eta$true_mosaicism * 100, "%")

  # Focus on the diagnostically relevant mosaicism range
  plot_data <- results_eta[results_eta$true_mosaicism %in% c(0.10, 0.15, 0.20), ]

  p <- ggplot(plot_data,
              aes(x = factor(coverage), y = prop_normal,
                  colour = eta_label, group = eta_label)) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, alpha = 0.5) +
    facet_wrap(~m_label) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    scale_colour_brewer(palette = "Dark2") +
    labs(
      x      = "Sequencing coverage (×)",
      y      = "Proportion classified as normal",
      colour = "PCR multiplier (η)",
      title  = "Supplementary: sensitivity to PCR noise multiplier (η)",
      subtitle = "Thresholds 0.15 / 0.30; n = 5,000 per cell"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave("figures/supp_eta_sensitivity.png", p,
         width = 9, height = 4, dpi = 300)
  cat("Saved: figures/supp_eta_sensitivity.png\n")
} else {
  cat("ggplot2 not available — skipping figure.\n")
}

# ── Console summary ─────────────────────────────────────────────────────────
cat("\n=== η sensitivity: false-normal rate at M = 15%, coverage = 10× ===\n")
sub <- results_eta[results_eta$true_mosaicism == 0.15 & results_eta$coverage == 10, ]
sub$range_pp <- round((sub$ci_upper - sub$ci_lower) * 100, 1)
print(sub[, c("eta", "prop_normal", "ci_lower", "ci_upper")],
      row.names = FALSE, digits = 3)

cat("\nMax absolute difference in prop_normal across η values (M=15%, 10×):",
    round(max(sub$prop_normal) - min(sub$prop_normal), 3), "\n")
