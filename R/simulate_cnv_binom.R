# simulate_cnv_binom.R
#
# Extension of simulate_cnv.R that adds binomial sampling of cells
# from the trophectoderm biopsy before computing the observed log2 ratio.
#
# Motivation: the main simulation treats true mosaicism M as fixed
# when computing the observed copy ratio. In reality, a biopsy of
# n_cells cells samples a random subset of the embryo: the number
# of cells carrying the CNV is Binomial(n_cells, M), so the
# *measured* mosaicism can deviate substantially from M by chance alone.
#
# Effect: at M = 0.15, 10x coverage, n_cells = 5, the false-normal rate
# rises from 71.5% (fixed-M model) to approximately 84% (binomial model).
# This is reported in the paper's Limitations section but was not included
# in the main analysis. This file makes it the base case.
#
# Author: Anere Cye

simulate_cnv_binom <- function(true_mosaicism,
                               cnv_type        = "dup",
                               coverage        = 10,
                               noise_sd_base   = 0.05,
                               pcr_duplicates  = 1.2,
                               n_cells         = 5) {
  # ── 1. Binomial sampling of biopsied cells ──────────────────────────────
  # Among n_cells biopsied, how many carry the CNV?
  # When n_cells is NA or 0, fall back to the fixed-M model (original behaviour)
  if (!is.na(n_cells) && n_cells > 0) {
    n_affected   <- rbinom(1, size = n_cells, prob = true_mosaicism)
    sampled_M    <- n_affected / n_cells
  } else {
    sampled_M    <- true_mosaicism
  }

  # ── 2. True copy ratio from sampled mosaicism ───────────────────────────
  if (cnv_type == "dup") {
    true_cr <- 1 + sampled_M / 2
  } else if (cnv_type == "del") {
    true_cr <- 1 - sampled_M / 2
  } else {
    stop("cnv_type must be 'dup' or 'del'")
  }

  # ── 3. Technical noise (same model as simulate_cnv.R) ───────────────────
  effective_noise_sd <- noise_sd_base * pcr_duplicates / sqrt(coverage / 10)
  observed_cr        <- true_cr + rnorm(1, mean = 0, sd = effective_noise_sd)

  # Guard against non-positive copy ratios (rare but possible at high noise)
  observed_cr <- max(observed_cr, 1e-6)

  log2_ratio <- log2(observed_cr)
  return(log2_ratio)
}
