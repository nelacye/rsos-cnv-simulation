# sensitivity_filter_factor.R
#
# Explicit sensitivity analysis for the post-calling filter reduction factor
# applied to false positive rates in Scenario B.
#
# In the paper (Section 2.3c and Table 4), the raw false positive rate from
# Scenario B (54.4%) is adjusted by a factor of 2–3 derived from published
# re-biopsy studies, yielding an estimated 18–27%. This factor is applied as
# a point estimate rather than as a range, and is not varied in the main code.
#
# This script:
#   1. Re-runs Scenario B to obtain the raw FP rate (with CI).
#   2. Applies filter factors of 1 (unfiltered), 2, 3, 5 to the raw rate.
#   3. Propagates uncertainty through the CI using delta method approximation.
#   4. Saves a table and an optional figure showing the trade-off between
#      FP reduction and the implied false-negative increase.
#
# Note on the circularity caveat (Section 3.5 / 4.7 of the paper):
#   The filter factors come from re-biopsy confirmation rates, which are
#   themselves threshold-dependent. This script documents the range of
#   filter-adjusted FP estimates rather than asserting any single value.
#
# Output:
#   results/filter_factor_sensitivity.csv
#   figures/supp_filter_factor.png   (if ggplot2 present)
#
# Author: Anere Cye

source("R/simulate_embryo_single_affected.R")
source("R/classify_embryo.R")

set.seed(2025)

# ── Parameters ─────────────────────────────────────────────────────────────
n_replicates   <- 5000
low_thr        <- 0.15
high_thr       <- 0.30
coverage       <- 10
filter_factors <- c(1, 2, 3, 5)   # 1 = no filtering

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ── Wilson CI helper ────────────────────────────────────────────────────────
wilson_ci <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  ci <- prop.test(k, n, conf.level = conf, correct = FALSE)$conf.int
  as.numeric(ci)
}

# ── Step 1: Raw Scenario B FP rate (M = 0) ─────────────────────────────────
cat("Running Scenario B (M = 0, n =", n_replicates, ")...\n")

fp_count <- 0L
for (i in seq_len(n_replicates)) {
  obs <- simulate_embryo_single_affected(0, coverage = coverage)
  cls <- classify_embryo(obs, low_thr, high_thr)
  if (any(cls != "normal")) fp_count <- fp_count + 1L
}

fp_raw <- fp_count / n_replicates
ci_raw <- wilson_ci(fp_count, n_replicates)

cat(sprintf("Raw FP rate: %.1f%% [%.1f%%–%.1f%%]\n",
            fp_raw * 100, ci_raw[1] * 100, ci_raw[2] * 100))

# ── Step 2: Apply filter factors and propagate uncertainty ──────────────────
# Delta-method approximation: if FP_adj = FP_raw / f, then
#   SE(FP_adj) ≈ SE(FP_raw) / f
# The Wilson CI is asymmetric; we apply the factor to both bounds
# as an approximation, clamping at [0, 1].

rows <- list()
for (f in filter_factors) {
  fp_adj     <- min(fp_raw / f, 1)
  ci_adj_lo  <- max(ci_raw[1] / f, 0)
  ci_adj_hi  <- min(ci_raw[2] / f, 1)

  rows[[length(rows) + 1]] <- data.frame(
    filter_factor       = f,
    fp_raw              = fp_raw,
    fp_adjusted         = fp_adj,
    ci_lower_adjusted   = ci_adj_lo,
    ci_upper_adjusted   = ci_adj_hi,
    filter_description  = switch(as.character(f),
      "1" = "No filter (upper bound)",
      "2" = "Conservative filter (lower end of published range)",
      "3" = "Moderate filter (upper end of published range)",
      "5" = "Aggressive filter"
    )
  )
}

results_ff <- do.call(rbind, rows)
write.csv(results_ff, "results/filter_factor_sensitivity.csv", row.names = FALSE)
cat("Saved: results/filter_factor_sensitivity.csv\n")

# ── Console table ───────────────────────────────────────────────────────────
cat("\n=== Filter-adjusted false positive rates (Scenario B, M = 0, 10×) ===\n")
cat(sprintf("%-20s  %-12s  %-20s\n", "Filter factor", "FP adj (%)", "95% CI (%)"))
cat(strrep("-", 58), "\n")
for (i in seq_len(nrow(results_ff))) {
  r <- results_ff[i, ]
  cat(sprintf("%-20s  %-12.1f  [%.1f%% – %.1f%%]\n",
              paste0("÷", r$filter_factor),
              r$fp_adjusted * 100,
              r$ci_lower_adjusted * 100,
              r$ci_upper_adjusted * 100))
}
cat("\nNote: CI propagation uses delta-method approximation on Wilson bounds.\n")
cat("Filter factors 2–3 correspond to published re-biopsy data\n")
cat("(Capalbo et al., 2021; Viotti et al., 2021).\n")

# ── Optional figure ─────────────────────────────────────────────────────────
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  results_ff$label <- factor(
    paste0("÷", results_ff$filter_factor),
    levels = paste0("÷", filter_factors)
  )

  # Shade the published range (factors 2–3)
  published_lo <- min(results_ff$fp_adjusted[results_ff$filter_factor %in% c(2, 3)])
  published_hi <- max(results_ff$fp_adjusted[results_ff$filter_factor %in% c(2, 3)])

  p <- ggplot(results_ff,
              aes(x = label, y = fp_adjusted * 100)) +
    annotate("rect",
             xmin = 1.5, xmax = 2.5,
             ymin = 0, ymax = 100,
             fill = "steelblue", alpha = 0.12) +
    annotate("text", x = 2, y = 95, size = 3, colour = "steelblue4",
             label = "Published\nfilter range") +
    geom_col(fill = "tomato3", width = 0.55, alpha = 0.85) +
    geom_errorbar(aes(ymin = ci_lower_adjusted * 100,
                      ymax = ci_upper_adjusted * 100),
                  width = 0.25) +
    geom_hline(yintercept = fp_raw * 100, linetype = "dashed", colour = "grey50") +
    annotate("text", x = 0.6, y = fp_raw * 100 + 2.5, size = 3,
             colour = "grey40", label = "Unfiltered") +
    scale_y_continuous(limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(
      x        = "Filter reduction factor",
      y        = "Filter-adjusted false positive rate",
      title    = "Supplementary: sensitivity of Scenario B FP rate to filter factor",
      subtitle = paste0("Raw FP = ",
                        round(fp_raw * 100, 1), "% [",
                        round(ci_raw[1] * 100, 1), "–",
                        round(ci_raw[2] * 100, 1), "%]; ",
                        "coverage = 10×; thresholds 0.15/0.30; n = 5,000")
    ) +
    theme_minimal(base_size = 11)

  ggsave("figures/supp_filter_factor.png", p,
         width = 7, height = 5, dpi = 300)
  cat("Saved: figures/supp_filter_factor.png\n")
} else {
  cat("ggplot2 not available — skipping figure.\n")
}
