# Autosomes-only validation on the pseudo log2-like workbook-derived signal.
#
# This script:
#   1. Filters to chr1-chr22
#   2. Defines an empirical chromosome abnormality cutoff from the autosomal
#      absolute log2-like distribution
#   3. Aggregates chromosome-level calls to embryo-level summaries
#   4. Compares embryo-level calls against EmbryoStatus

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

run_wgspgt_autosomal_validation <- function(
    input_path = file.path("results", "wgspgt_descriptive", "wgspgt_log2like.csv"),
    results_dir = file.path("results", "wgspgt_autosomal_validation"),
    quantile_cutoff = 0.95) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- readr::read_csv(input_path, show_col_types = FALSE) %>%
    mutate(
      embryo_id = as.character(embryo_id),
      chromosome = as.character(chromosome),
      EmbryoStatus = as.character(EmbryoStatus),
      Method = as.character(Method)
    )

  autosome_levels <- paste0("chr", 1:22)

  autosomal <- dat %>%
    filter(chromosome %in% autosome_levels)

  chromosome_cutoff <- quantile(abs(autosomal$log2like), probs = quantile_cutoff, na.rm = TRUE)

  chromosome_calls <- autosomal %>%
    mutate(
      abs_log2like = abs(log2like),
      chromosome_call = if_else(abs_log2like > chromosome_cutoff, "abnormal", "normal")
    )

  embryo_level <- chromosome_calls %>%
    group_by(embryo_id, EmbryoStatus, Method, meanDepth, meanBreadth) %>%
    summarise(
      n_autosomes = n(),
      n_abnormal_chromosomes = sum(chromosome_call == "abnormal", na.rm = TRUE),
      mean_log2like = mean(log2like, na.rm = TRUE),
      variance_log2like = var(log2like, na.rm = TRUE),
      max_abs_log2like = max(abs_log2like, na.rm = TRUE),
      embryo_call = if_else(n_abnormal_chromosomes > 0, "abnormal", "normal"),
      .groups = "drop"
    ) %>%
    mutate(
      truth_binary = if_else(EmbryoStatus == "unaffected", "normal", "abnormal")
    ) %>%
    arrange(as.numeric(embryo_id))

  confusion_matrix <- embryo_level %>%
    count(truth_binary, embryo_call, name = "n") %>%
    arrange(truth_binary, embryo_call)

  tp <- sum(embryo_level$truth_binary == "abnormal" & embryo_level$embryo_call == "abnormal")
  tn <- sum(embryo_level$truth_binary == "normal" & embryo_level$embryo_call == "normal")
  fp <- sum(embryo_level$truth_binary == "normal" & embryo_level$embryo_call == "abnormal")
  fn <- sum(embryo_level$truth_binary == "abnormal" & embryo_level$embryo_call == "normal")

  metrics <- tibble(
    quantile_cutoff = quantile_cutoff,
    chromosome_abs_cutoff = as.numeric(chromosome_cutoff),
    n_embryos = nrow(embryo_level),
    sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
    specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
    accuracy = (tp + tn) / nrow(embryo_level),
    true_positive = tp,
    true_negative = tn,
    false_positive = fp,
    false_negative = fn
  )

  by_status <- embryo_level %>%
    group_by(EmbryoStatus, embryo_call) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(EmbryoStatus, embryo_call)

  readr::write_csv(chromosome_calls, file.path(results_dir, "autosomal_chromosome_calls.csv"))
  readr::write_csv(embryo_level, file.path(results_dir, "per_embryo_classification.csv"))
  readr::write_csv(confusion_matrix, file.path(results_dir, "confusion_matrix.csv"))
  readr::write_csv(metrics, file.path(results_dir, "validation_metrics.csv"))
  readr::write_csv(by_status, file.path(results_dir, "classification_by_embryo_status.csv"))

  report_lines <- c(
    "# Autosomal Validation Report",
    "",
    "## Scope",
    "",
    "This analysis uses `wgspgt_log2like.csv` and restricts validation to autosomes only (`chr1`-`chr22`).",
    "",
    "## Empirical Threshold",
    "",
    sprintf("- Chromosome abnormality threshold: `|log2like| > %.3f`", chromosome_cutoff),
    sprintf("- Threshold source: %.0fth percentile of `|log2like|` across autosomes", quantile_cutoff * 100),
    "",
    "## Embryo-Level Rule",
    "",
    "- An embryo is called `abnormal` if it has at least one abnormal autosome.",
    "- Embryo-level summaries include number of abnormal autosomes, mean log2-like value, and variance.",
    "",
    "## EmbryoStatus Comparison",
    "",
    "- Binary truth mapping used for metrics:",
    "  - `unaffected` -> `normal`",
    "  - all other statuses -> `abnormal`",
    "- This is a pragmatic first-pass comparison because no dedicated diagnosis sheet is being used here.",
    "",
    "## Performance",
    "",
    sprintf("- Sensitivity: %.3f", metrics$sensitivity[[1]]),
    sprintf("- Specificity: %.3f", metrics$specificity[[1]]),
    sprintf("- Accuracy: %.3f", metrics$accuracy[[1]]),
    sprintf("- TP / TN / FP / FN: %s / %s / %s / %s", tp, tn, fp, fn),
    "",
    "## Interpretation",
    "",
    "- These results should be treated as exploratory because the pseudo log2-like transformation is not yet equivalent to original QDNAseq logR.",
    "- The binary truth mapping based on `EmbryoStatus` is also only an approximation of CNV abnormality.",
    "- Even so, this autosomes-only pass is a useful checkpoint for whether the transformed signal carries embryo-level abnormality information."
  )

  writeLines(report_lines, file.path(results_dir, "autosomal_validation_report.md"))

  message("WGSPGT autosomal validation complete.")
  message("  Results directory: ", results_dir)
  message("  Per-embryo table: ", file.path(results_dir, "per_embryo_classification.csv"))
  message("  Report: ", file.path(results_dir, "autosomal_validation_report.md"))

  invisible(
    list(
      chromosome_calls = chromosome_calls,
      embryo_level = embryo_level,
      confusion_matrix = confusion_matrix,
      metrics = metrics
    )
  )
}

if (sys.nframe() == 0L) {
  run_wgspgt_autosomal_validation()
}
