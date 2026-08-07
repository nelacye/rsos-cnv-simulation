# Compare autosomal validation performance under strict and inclusive
# EmbryoStatus mappings.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

run_wgspgt_status_mapping_validation <- function(
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

  mapping_table <- tibble(
    EmbryoStatus = c("affected", "unaffected", "carrier", "inconclusive"),
    clinical_group = c("clearly_abnormal", "clearly_normal", "ambiguous_intermediate", "ambiguous_intermediate"),
    strict_truth = c("abnormal", "normal", NA, NA),
    inclusive_truth = c("abnormal", "normal", "abnormal", "abnormal")
  )

  readr::write_csv(mapping_table, file.path(results_dir, "embryo_status_mapping.csv"))

  autosomal <- dat %>%
    filter(chromosome %in% autosome_levels) %>%
    left_join(mapping_table, by = "EmbryoStatus")

  chromosome_cutoff <- quantile(abs(autosomal$log2like), probs = quantile_cutoff, na.rm = TRUE)

  chromosome_calls <- autosomal %>%
    mutate(
      abs_log2like = abs(log2like),
      chromosome_call = if_else(abs_log2like > chromosome_cutoff, "abnormal", "normal")
    )

  embryo_level <- chromosome_calls %>%
    group_by(embryo_id, EmbryoStatus, clinical_group, strict_truth, inclusive_truth, Method, meanDepth, meanBreadth) %>%
    summarise(
      n_autosomes = n(),
      n_abnormal_chromosomes = sum(chromosome_call == "abnormal", na.rm = TRUE),
      mean_log2like = mean(log2like, na.rm = TRUE),
      variance_log2like = var(log2like, na.rm = TRUE),
      max_abs_log2like = max(abs_log2like, na.rm = TRUE),
      embryo_call = if_else(n_abnormal_chromosomes > 0, "abnormal", "normal"),
      .groups = "drop"
    ) %>%
    arrange(as.numeric(embryo_id))

  evaluate_mode <- function(df, truth_col, mode_name) {
    eval_df <- df %>%
      filter(!is.na(.data[[truth_col]])) %>%
      mutate(truth_binary = .data[[truth_col]])

    tp <- sum(eval_df$truth_binary == "abnormal" & eval_df$embryo_call == "abnormal")
    tn <- sum(eval_df$truth_binary == "normal" & eval_df$embryo_call == "normal")
    fp <- sum(eval_df$truth_binary == "normal" & eval_df$embryo_call == "abnormal")
    fn <- sum(eval_df$truth_binary == "abnormal" & eval_df$embryo_call == "normal")

    list(
      confusion = eval_df %>%
        count(truth_binary, embryo_call, name = "n") %>%
        mutate(mode = mode_name) %>%
        select(mode, truth_binary, embryo_call, n) %>%
        arrange(mode, truth_binary, embryo_call),
      metrics = tibble(
        mode = mode_name,
        quantile_cutoff = quantile_cutoff,
        chromosome_abs_cutoff = as.numeric(chromosome_cutoff),
        n_embryos = nrow(eval_df),
        sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
        specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
        accuracy = (tp + tn) / nrow(eval_df),
        true_positive = tp,
        true_negative = tn,
        false_positive = fp,
        false_negative = fn
      ),
      per_embryo = eval_df %>%
        mutate(mode = mode_name) %>%
        select(mode, embryo_id, EmbryoStatus, clinical_group, truth_binary, embryo_call,
               n_autosomes, n_abnormal_chromosomes, mean_log2like, variance_log2like,
               max_abs_log2like, Method, meanDepth, meanBreadth)
    )
  }

  strict_eval <- evaluate_mode(embryo_level, "strict_truth", "strict")
  inclusive_eval <- evaluate_mode(embryo_level, "inclusive_truth", "inclusive")

  metrics_comparison <- bind_rows(strict_eval$metrics, inclusive_eval$metrics)
  confusion_comparison <- bind_rows(strict_eval$confusion, inclusive_eval$confusion)
  per_embryo_comparison <- bind_rows(strict_eval$per_embryo, inclusive_eval$per_embryo)

  readr::write_csv(confusion_comparison, file.path(results_dir, "confusion_matrix_by_mode.csv"))
  readr::write_csv(metrics_comparison, file.path(results_dir, "validation_metrics_by_mode.csv"))
  readr::write_csv(per_embryo_comparison, file.path(results_dir, "per_embryo_classification_by_mode.csv"))

  report_lines <- c(
    "# EmbryoStatus Mapping Validation Report",
    "",
    "## Unique EmbryoStatus Values",
    "",
    "- affected",
    "- carrier",
    "- inconclusive",
    "- unaffected",
    "",
    "## Clinical Grouping",
    "",
    "- clearly abnormal: `affected`",
    "- clearly normal: `unaffected`",
    "- ambiguous / intermediate: `carrier`, `inconclusive`",
    "",
    "## Validation Modes",
    "",
    sprintf("- Shared chromosome abnormality threshold: `|log2like| > %.3f`", chromosome_cutoff),
    sprintf("- Threshold source: %.0fth percentile of autosomal absolute log2-like values", quantile_cutoff * 100),
    "",
    "### Strict",
    "",
    "- Includes only clear normal vs clear abnormal embryos.",
    "- Excludes ambiguous/intermediate statuses.",
    sprintf("- Sensitivity: %.3f", strict_eval$metrics$sensitivity[[1]]),
    sprintf("- Specificity: %.3f", strict_eval$metrics$specificity[[1]]),
    sprintf("- Accuracy: %.3f", strict_eval$metrics$accuracy[[1]]),
    sprintf("- TP / TN / FP / FN: %s / %s / %s / %s",
            strict_eval$metrics$true_positive[[1]],
            strict_eval$metrics$true_negative[[1]],
            strict_eval$metrics$false_positive[[1]],
            strict_eval$metrics$false_negative[[1]]),
    "",
    "### Inclusive",
    "",
    "- Uses the current mapping where non-`unaffected` statuses are treated as abnormal.",
    sprintf("- Sensitivity: %.3f", inclusive_eval$metrics$sensitivity[[1]]),
    sprintf("- Specificity: %.3f", inclusive_eval$metrics$specificity[[1]]),
    sprintf("- Accuracy: %.3f", inclusive_eval$metrics$accuracy[[1]]),
    sprintf("- TP / TN / FP / FN: %s / %s / %s / %s",
            inclusive_eval$metrics$true_positive[[1]],
            inclusive_eval$metrics$true_negative[[1]],
            inclusive_eval$metrics$false_positive[[1]],
            inclusive_eval$metrics$false_negative[[1]]),
    "",
    "## Comparison",
    "",
    "- The strict mode is clinically cleaner because it removes ambiguous statuses from the truth set.",
    "- The inclusive mode keeps more embryos but mixes intermediate categories into the abnormal class.",
    "- Use `validation_metrics_by_mode.csv` and `per_embryo_classification_by_mode.csv` for side-by-side comparison."
  )

  writeLines(report_lines, file.path(results_dir, "status_mapping_validation_report.md"))

  message("WGSPGT status-mapping validation comparison complete.")
  message("  Mapping table: ", file.path(results_dir, "embryo_status_mapping.csv"))
  message("  Metrics: ", file.path(results_dir, "validation_metrics_by_mode.csv"))
  message("  Report: ", file.path(results_dir, "status_mapping_validation_report.md"))

  invisible(
    list(
      mapping_table = mapping_table,
      metrics_comparison = metrics_comparison,
      confusion_comparison = confusion_comparison,
      per_embryo_comparison = per_embryo_comparison
    )
  )
}

if (sys.nframe() == 0L) {
  run_wgspgt_status_mapping_validation()
}
