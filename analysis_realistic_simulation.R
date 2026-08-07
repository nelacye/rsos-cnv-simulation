# Prototype realistic aggregated-signal simulation.
#
# This script keeps the original idealized simulation untouched and adds a
# parallel realistic workflow that models:
#   1. latent embryo state
#   2. latent chromosome signal
#   3. measurement distortion / aggregation
#   4. embryo-level decision rules

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

source("R/simulate_latent_embryo_state.R")
source("R/simulate_latent_chromosome_signal.R")
source("R/simulate_measurement_layer.R")
source("R/classify_embryo_realistic.R")
source("R/simulate_embryo_single_affected.R")
source("R/simulate_embryo_comparable_idealized.R")
source("R/classify_embryo.R")

run_idealized_comparison <- function(latent_states, low_thr = 0.15, high_thr = 0.30) {
  status_to_mosaicism <- c(normal = 0.00, ambiguous = 0.15, abnormal = 0.25)

  embryo_eval <- vector("list", nrow(latent_states))

  for (i in seq_len(nrow(latent_states))) {
    status_i <- latent_states$latent_status[i]
    m_i <- status_to_mosaicism[[status_i]]
    obs_i <- simulate_embryo_single_affected(
      true_mosaicism = m_i,
      n_chr = 22,
      coverage = latent_states$coverage[i]
    )
    classes_i <- classify_embryo(obs_i, low_thr = low_thr, high_thr = high_thr)
    embryo_call_i <- ifelse(any(classes_i != "normal"), "abnormal", "normal")
    truth_i <- ifelse(status_i == "normal", "normal", "abnormal")

    embryo_eval[[i]] <- data.frame(
      embryo_id = latent_states$embryo_id[i],
      latent_status = status_i,
      classifier = "idealized_threshold",
      rule_name = "ge_1_abnormal_chromosome",
      embryo_call = embryo_call_i,
      truth_binary = truth_i,
      stringsAsFactors = FALSE
    )
  }

  idealized_calls <- do.call(rbind, embryo_eval)

  tp <- sum(idealized_calls$truth_binary == "abnormal" & idealized_calls$embryo_call == "abnormal")
  tn <- sum(idealized_calls$truth_binary == "normal" & idealized_calls$embryo_call == "normal")
  fp <- sum(idealized_calls$truth_binary == "normal" & idealized_calls$embryo_call == "abnormal")
  fn <- sum(idealized_calls$truth_binary == "abnormal" & idealized_calls$embryo_call == "normal")

  metrics <- data.frame(
    classifier = "idealized_threshold",
    rule_name = "ge_1_abnormal_chromosome",
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    accuracy = (tp + tn) / nrow(idealized_calls),
    true_positive = tp,
    true_negative = tn,
    false_positive = fp,
    false_negative = fn,
    stringsAsFactors = FALSE
  )

  list(
    calls = idealized_calls,
    metrics = metrics
  )
}

run_comparable_idealized_comparison <- function(latent_states, low_thr = 0.15, high_thr = 0.30) {
  status_to_mosaicism <- c(normal = 0.00, ambiguous = 0.15, abnormal = 0.25)
  status_to_n_affected <- list(
    normal = 0L,
    ambiguous = c(1L, 2L),
    abnormal = 2:6
  )

  embryo_eval <- vector("list", nrow(latent_states))

  for (i in seq_len(nrow(latent_states))) {
    status_i <- latent_states$latent_status[i]
    m_i <- status_to_mosaicism[[status_i]]
    n_affected_i <- if (length(status_to_n_affected[[status_i]]) == 1L) {
      status_to_n_affected[[status_i]]
    } else {
      sample(status_to_n_affected[[status_i]], size = 1)
    }

    obs_i <- simulate_embryo_comparable_idealized(
      true_mosaicism = m_i,
      n_affected = n_affected_i,
      n_chr = 22,
      coverage = latent_states$coverage[i]
    )
    classes_i <- classify_embryo(obs_i, low_thr = low_thr, high_thr = high_thr)
    embryo_call_i <- ifelse(any(classes_i != "normal"), "abnormal", "normal")
    truth_i <- ifelse(status_i == "normal", "normal", "abnormal")

    embryo_eval[[i]] <- data.frame(
      embryo_id = latent_states$embryo_id[i],
      latent_status = status_i,
      n_affected = n_affected_i,
      classifier = "comparable_idealized_threshold",
      rule_name = "ge_1_abnormal_chromosome",
      embryo_call = embryo_call_i,
      truth_binary = truth_i,
      stringsAsFactors = FALSE
    )
  }

  idealized_calls <- do.call(rbind, embryo_eval)

  tp <- sum(idealized_calls$truth_binary == "abnormal" & idealized_calls$embryo_call == "abnormal")
  tn <- sum(idealized_calls$truth_binary == "normal" & idealized_calls$embryo_call == "normal")
  fp <- sum(idealized_calls$truth_binary == "normal" & idealized_calls$embryo_call == "abnormal")
  fn <- sum(idealized_calls$truth_binary == "abnormal" & idealized_calls$embryo_call == "normal")

  metrics <- data.frame(
    classifier = "comparable_idealized_threshold",
    rule_name = "ge_1_abnormal_chromosome",
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    accuracy = (tp + tn) / nrow(idealized_calls),
    true_positive = tp,
    true_negative = tn,
    false_positive = fp,
    false_negative = fn,
    stringsAsFactors = FALSE
  )

  list(
    calls = idealized_calls,
    metrics = metrics
  )
}

run_realistic_simulation <- function(
    n_embryos = 1500,
    seed = 2026,
    calibration_fraction = 0.50,
    results_dir = file.path("results", "realistic_simulation")) {
  set.seed(seed)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  latent_states <- simulate_latent_embryo_state(n_embryos = n_embryos)
  latent_chromosomes <- simulate_latent_chromosome_signal(latent_states)
  observed_chromosomes <- simulate_measurement_layer(
    latent_chromosome_signal = latent_chromosomes,
    embryo_state = latent_states
  )

  chromosome_cutoff <- estimate_empirical_signal_cutoff(observed_chromosomes, quantile_cutoff = 0.95)
  chromosome_calls <- compute_chromosome_calls_realistic(
    observed_data = observed_chromosomes,
    chromosome_cutoff = chromosome_cutoff,
    posterior_cutoff = 0.90
  )

  embryo_summary <- summarize_embryo_realistic(chromosome_calls)
  embryo_ids <- unique(embryo_summary$embryo_id)
  n_calibration <- floor(length(embryo_ids) * calibration_fraction)

  calibration_ids <- sample(embryo_ids, size = n_calibration, replace = FALSE)
  embryo_summary$data_split <- ifelse(
    embryo_summary$embryo_id %in% calibration_ids,
    "calibration",
    "evaluation"
  )

  calibration_summary <- embryo_summary[embryo_summary$data_split == "calibration", ]
  evaluation_summary <- embryo_summary[embryo_summary$data_split == "evaluation", ]

  calibration_cutoffs <- estimate_aggregate_cutoffs_realistic(
    calibration_summary = calibration_summary,
    evaluation_label = "out-of-sample evaluation"
  )

  decision_results <- apply_embryo_decision_rules_realistic(
    embryo_summary = evaluation_summary,
    aggregate_cutoffs = calibration_cutoffs,
    evaluation_label = "out-of-sample evaluation"
  )
  scored <- score_realistic_decisions(decision_results$embryo_calls)
  idealized <- run_idealized_comparison(latent_states)
  comparable_idealized <- run_comparable_idealized_comparison(latent_states)

  realistic_vs_idealized <- rbind(
    scored$metrics[, c("classifier", "rule_name", "sensitivity", "specificity", "accuracy",
                       "true_positive", "true_negative", "false_positive", "false_negative")],
    comparable_idealized$metrics
  )

  signal_hist <- ggplot(chromosome_calls, aes(x = observed_signal, fill = latent_status)) +
    geom_histogram(position = "identity", alpha = 0.35, bins = 60) +
    theme_minimal(base_size = 11) +
    labs(
      title = "Realistic observed chromosome-level signal",
      x = "Observed pseudo signal",
      y = "Count",
      fill = "Latent status"
    )

  signal_density <- ggplot(chromosome_calls, aes(x = observed_signal, color = latent_status)) +
    geom_density(linewidth = 0.8) +
    theme_minimal(base_size = 11) +
    labs(
      title = "Observed chromosome-level signal by latent embryo state",
      x = "Observed pseudo signal",
      y = "Density",
      color = "Latent status"
    )

  comparison_plot <- ggplot(realistic_vs_idealized, aes(x = rule_name, y = specificity, fill = classifier)) +
    geom_col(position = "dodge") +
    geom_point(aes(y = sensitivity, color = classifier), position = position_dodge(width = 0.9), size = 2) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(
      title = "Idealized vs realistic simulation outputs",
      subtitle = "Realistic results use out-of-sample evaluation; bars show specificity, points show sensitivity",
      x = "Decision rule",
      y = "Performance",
      fill = "Classifier",
      color = "Classifier"
    )

  ggsave(file.path(results_dir, "observed_signal_histogram.png"), signal_hist, width = 8, height = 5, dpi = 300)
  ggsave(file.path(results_dir, "observed_signal_density.png"), signal_density, width = 8, height = 5, dpi = 300)
  ggsave(file.path(results_dir, "idealized_vs_realistic_comparison.png"), comparison_plot, width = 10, height = 5, dpi = 300)

  write_csv(latent_states, file.path(results_dir, "latent_embryo_states.csv"))
  write_csv(latent_chromosomes, file.path(results_dir, "latent_chromosome_signal.csv"))
  write_csv(chromosome_calls, file.path(results_dir, "observed_chromosome_signal.csv"))
  write_csv(embryo_summary, file.path(results_dir, "embryo_signal_summary.csv"))
  write_csv(calibration_cutoffs, file.path(results_dir, "calibration_cutoffs.csv"))
  write_csv(calibration_cutoffs, file.path(results_dir, "aggregate_rule_thresholds.csv"))
  write_csv(decision_results$embryo_calls, file.path(results_dir, "embryo_decision_calls.csv"))
  write_csv(scored$metrics, file.path(results_dir, "evaluation_metrics.csv"))
  write_csv(scored$metrics, file.path(results_dir, "decision_rule_comparison.csv"))
  write_csv(scored$confusion, file.path(results_dir, "confusion_matrices.csv"))
  write_csv(realistic_vs_idealized, file.path(results_dir, "idealized_vs_realistic_metrics.csv"))
  write_csv(idealized$calls, file.path(results_dir, "idealized_embryo_calls.csv"))
  write_csv(comparable_idealized$calls, file.path(results_dir, "comparable_idealized_embryo_calls.csv"))

  report_lines <- c(
    "# Realistic Simulation Report",
    "",
    "## Scope",
    "",
    sprintf("- Simulated embryos: %s", n_embryos),
    sprintf("- Calibration fraction: %.2f", calibration_fraction),
    sprintf("- Empirical chromosome cutoff for fixed-threshold logic: %.3f", chromosome_cutoff),
    "- Latent statuses: normal, ambiguous, abnormal",
    "- Autosomes only (`chr1`-`chr22`)",
    "- Aggregate embryo rules are reported as out-of-sample evaluation",
    "",
    "## Interpretation",
    "",
    "- The realistic model inserts an explicit measurement layer between latent chromosome state and observed pseudo-signal.",
    "- Aggregate embryo cutoffs are estimated on a calibration split using only latent normal embryos.",
    "- Threshold and Bayesian chromosome logic are evaluated on the held-out evaluation split.",
    "- The idealized-vs-realistic comparison now uses a parallel comparable-idealized mode with matched latent chromosome burden, so it isolates the measurement layer more cleanly.",
    "- The output tables can be compared directly through sensitivity, specificity, accuracy, and confusion counts.",
    "",
    "## Key Output Files",
    "",
    "- `decision_rule_comparison.csv`",
    "- `evaluation_metrics.csv`",
    "- `calibration_cutoffs.csv`",
    "- `confusion_matrices.csv`",
    "- `idealized_vs_realistic_metrics.csv`",
    "- `observed_chromosome_signal.csv`",
    "- `embryo_decision_calls.csv`"
  )

  writeLines(report_lines, file.path(results_dir, "realistic_simulation_report.md"))

  invisible(
    list(
      latent_states = latent_states,
      chromosome_calls = chromosome_calls,
      embryo_summary = embryo_summary,
      calibration_cutoffs = calibration_cutoffs,
      decision_metrics = scored$metrics,
      idealized_metrics = idealized$metrics,
      comparable_idealized_metrics = comparable_idealized$metrics
    )
  )
}

if (sys.nframe() == 0L) {
  run_realistic_simulation()
}
