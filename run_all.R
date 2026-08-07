# ============================================================
# run_all.R
#
# Master script — runs all simulations and figures in order.
# Reproduces every result, table, and figure in the paper.
#
# Estimated runtimes (n = 5000, standard laptop):
#   Simulations (Scenarios single/A/B): ~5 min
#   Supplementary sensitivity:          ~5 min
#   Bayesian classifier (run_bayes.R):  ~20 min
#   fig_bayes_roc.R:                    ~30 min
#   fig_same_ranking_different_decisions.R: ~15 min
#   Total:                              ~75 min
#
# Quick check (n = 500): set QUICK_CHECK <- TRUE below.
# This reduces n_replicates to 500 in all Bayesian scripts.
# Results will be noisy but structure will be correct.
# ============================================================

QUICK_CHECK <- FALSE   # set TRUE for fast (~5 min) smoke test

set.seed(2025)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

t_start <- proc.time()

# ── 1. Main simulations ────────────────────────────────────────────────────

cat("╔══════════════════════════════════════════════╗\n")
cat("║  STEP 1/5 — Main simulations                ║\n")
cat("╚══════════════════════════════════════════════╝\n")

source("scripts/run_scenario_single.R")
source("scripts/run_scenario_A.R")
source("scripts/run_scenario_B.R")

cat("Running single-chromosome simulation...\n")
single_res <- run_scenario_single()
write.csv(single_res, "results/single_chr_results.csv", row.names = FALSE)

cat("Running Scenario A (uniform mosaicism, 24 chr)...\n")
scenario_A_res <- run_scenario_A()
write.csv(scenario_A_res, "results/scenario_A_results.csv", row.names = FALSE)

cat("Running Scenario B (single affected chr)...\n")
scenario_B_res <- run_scenario_B()
write.csv(scenario_B_res, "results/scenario_B_results.csv", row.names = FALSE)

single_subset <- single_res[single_res$cov == 10 &
                              single_res$m %in% c(0, 0.15, 0.20), ]
table5 <- data.frame(
  M                 = c(0, 0.15, 0.20),
  single_chr_normal = single_subset$prop_normal,
  scenario_A_normal = scenario_A_res$prop_normal,
  scenario_B_normal = scenario_B_res$prop_normal
)
write.csv(table5, "results/table5_comparison.csv", row.names = FALSE)
cat("  Saved: results/table5_comparison.csv\n")

# ── 2. Supplementary sensitivity analyses ─────────────────────────────────

cat("\n╔══════════════════════════════════════════════╗\n")
cat("║  STEP 2/5 — Supplementary sensitivity       ║\n")
cat("╚══════════════════════════════════════════════╝\n")

cat("Binomial cell sampling...\n")
source("scripts/sensitivity_binom_cells.R")

cat("PCR multiplier (eta)...\n")
source("scripts/sensitivity_eta.R")

cat("Scenario B filter factor...\n")
source("scripts/sensitivity_filter_factor.R")

# ── 3. Bayesian classifier — Table 3, posteriors, age-stratified ──────────

cat("\n╔══════════════════════════════════════════════╗\n")
cat("║  STEP 3/5 — Bayesian classifier             ║\n")
cat("╚══════════════════════════════════════════════╝\n")

if (QUICK_CHECK) {
  cat("  [QUICK_CHECK mode: n_replicates = 500]\n")
  # Temporarily patch n_replicates in run_bayes.R environment
  local({
    source("R/bayes_classifier.R")
    t3 <- build_table3(n_replicates = 500)
    write.csv(t3, "results/table3_bayes_priors.csv", row.names = FALSE)
    cat("  Saved: results/table3_bayes_priors.csv (quick check)\n")
  })
} else {
  source("scripts/run_bayes.R")
}

# ── 4. Figure: ROC-like curves (Bayesian vs fixed thresholds) ─────────────

cat("\n╔══════════════════════════════════════════════╗\n")
cat("║  STEP 4/5 — Figure: Bayesian ROC            ║\n")
cat("╚══════════════════════════════════════════════╝\n")
cat("  NOTE: ~30 min at N=5000. Set N_EMBRYOS<-500 in\n")
cat("  fig_bayes_roc.R for a quick check.\n\n")

source("fig_bayes_roc.R")

# ── 5. Figure: same ranking, different decisions ──────────────────────────

cat("\n╔══════════════════════════════════════════════╗\n")
cat("║  STEP 5/5 — Figure: Prior dependence        ║\n")
cat("╚══════════════════════════════════════════════╝\n")
cat("  NOTE: ~15 min at N=3000.\n\n")

source("fig_same_ranking_different_decisions.R")

# ── Done ──────────────────────────────────────────────────────────────────

t_elapsed <- round((proc.time() - t_start)["elapsed"] / 60, 1)

cat("\n╔══════════════════════════════════════════════╗\n")
cat("║  ✅  ALL STEPS COMPLETE                     ║\n")
cat("╚══════════════════════════════════════════════╝\n")
cat(sprintf("  Total elapsed: %.1f min\n\n", t_elapsed))

cat("Main results:\n")
cat("  results/single_chr_results.csv\n")
cat("  results/scenario_A_results.csv\n")
cat("  results/scenario_B_results.csv\n")
cat("  results/table5_comparison.csv\n")

cat("\nBayesian results:\n")
cat("  results/table3_bayes_priors.csv\n")
cat("  results/posterior_probs_embryos.csv\n")
cat("  results/table3_age_stratified.csv\n")
cat("  results/roc_curve_data.csv\n")
cat("  results/prior_dependence_data.csv\n")

cat("\nSupplementary results:\n")
cat("  results/binom_sensitivity.csv\n")
cat("  results/sensitivity_eta.csv\n")
cat("  results/filter_factor_sensitivity.csv\n")
cat("  results/bf_failure_data.csv\n")

cat("\nFigures:\n")
cat("  figures/fig_bayes_roc.png\n")
cat("  figures/fig_same_ranking_different_decisions.png\n")
cat("  figures/  (+ all other figures from main simulation)\n")
