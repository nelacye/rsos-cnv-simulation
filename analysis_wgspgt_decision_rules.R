# Compare alternative embryo-level decision rules on the autosomal
# pseudo log2-like signal under strict mode only.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

run_wgspgt_decision_rules <- function(
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

  status_mapping <- tibble(
    EmbryoStatus = c("affected", "unaffected", "carrier", "inconclusive"),
    strict_truth = c("abnormal", "normal", NA, NA)
  )

  autosomal <- dat %>%
    filter(chromosome %in% autosome_levels) %>%
    left_join(status_mapping, by = "EmbryoStatus")

  chromosome_abs_cutoff <- quantile(abs(autosomal$log2like), probs = quantile_cutoff, na.rm = TRUE)

  embryo_summary <- autosomal %>%
    group_by(embryo_id, EmbryoStatus, strict_truth, Method, meanDepth, meanBreadth) %>%
    summarise(
      n_autosomes = n(),
      n_abnormal_chromosomes = sum(abs(log2like) > chromosome_abs_cutoff, na.rm = TRUE),
      mean_abs_log2like = mean(abs(log2like), na.rm = TRUE),
      sum_abs_log2like = sum(abs(log2like), na.rm = TRUE),
      .groups = "drop"
    )

  strict_df <- embryo_summary %>%
    filter(!is.na(strict_truth))

  normal_reference <- strict_df %>%
    filter(strict_truth == "normal")

  mean_abs_threshold <- quantile(normal_reference$mean_abs_log2like, probs = 0.95, na.rm = TRUE)
  sum_abs_threshold <- quantile(normal_reference$sum_abs_log2like, probs = 0.95, na.rm = TRUE)

  rules <- list(
    list(name = "ge_1_abnormal_chromosome", type = "count", threshold = 1),
    list(name = "ge_2_abnormal_chromosomes", type = "count", threshold = 2),
    list(name = "ge_3_abnormal_chromosomes", type = "count", threshold = 3),
    list(name = "mean_abs_log2like_gt_threshold", type = "mean_abs", threshold = as.numeric(mean_abs_threshold)),
    list(name = "sum_abs_log2like_gt_threshold", type = "sum_abs", threshold = as.numeric(sum_abs_threshold))
  )

  evaluate_rule <- function(df, rule) {
    called <- if (rule$type == "count") {
      df %>%
        mutate(embryo_call = if_else(n_abnormal_chromosomes >= rule$threshold, "abnormal", "normal"))
    } else if (rule$type == "mean_abs") {
      df %>%
        mutate(embryo_call = if_else(mean_abs_log2like > rule$threshold, "abnormal", "normal"))
    } else if (rule$type == "sum_abs") {
      df %>%
        mutate(embryo_call = if_else(sum_abs_log2like > rule$threshold, "abnormal", "normal"))
    } else {
      stop(sprintf("Unsupported rule type: %s", rule$type), call. = FALSE)
    }

    tp <- sum(called$strict_truth == "abnormal" & called$embryo_call == "abnormal")
    tn <- sum(called$strict_truth == "normal" & called$embryo_call == "normal")
    fp <- sum(called$strict_truth == "normal" & called$embryo_call == "abnormal")
    fn <- sum(called$strict_truth == "abnormal" & called$embryo_call == "normal")

    confusion <- called %>%
      count(strict_truth, embryo_call, name = "n") %>%
      mutate(
        rule_name = rule$name,
        rule_threshold = rule$threshold
      ) %>%
      select(rule_name, rule_threshold, strict_truth, embryo_call, n) %>%
      arrange(rule_name, strict_truth, embryo_call)

    metrics <- tibble(
      rule_name = rule$name,
      rule_type = rule$type,
      rule_threshold = rule$threshold,
      chromosome_abs_cutoff = as.numeric(chromosome_abs_cutoff),
      n_embryos = nrow(called),
      sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
      specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
      accuracy = (tp + tn) / nrow(called),
      true_positive = tp,
      true_negative = tn,
      false_positive = fp,
      false_negative = fn
    )

    list(confusion = confusion, metrics = metrics)
  }

  evaluations <- lapply(rules, function(rule) evaluate_rule(strict_df, rule))

  confusion_all <- bind_rows(lapply(evaluations, `[[`, "confusion"))
  metrics_all <- bind_rows(lapply(evaluations, `[[`, "metrics")) %>%
    arrange(desc(specificity), desc(sensitivity), desc(accuracy))

  baseline <- metrics_all %>% filter(rule_name == "ge_1_abnormal_chromosome")
  improved_rules <- metrics_all %>%
    filter(
      specificity > baseline$specificity[[1]],
      sensitivity >= baseline$sensitivity[[1]] - 0.10
    )

  write_csv(confusion_all, file.path(results_dir, "decision_rule_confusion_matrices.csv"))
  write_csv(metrics_all, file.path(results_dir, "decision_rule_comparison.csv"))

  report_lines <- c(
    "# Decision Rule Comparison Report",
    "",
    "## Scope",
    "",
    "This comparison uses `wgspgt_log2like.csv`, autosomes only, and strict mode only.",
    "",
    "## Shared Inputs",
    "",
    sprintf("- Chromosome abnormality cutoff for count-based rules: `|log2like| > %.3f`", chromosome_abs_cutoff),
    sprintf("- Mean absolute log2like threshold: `%.3f`", mean_abs_threshold),
    sprintf("- Sum absolute log2like threshold: `%.3f`", sum_abs_threshold),
    "- Mean and sum thresholds were set from the 95th percentile of the clear-normal embryos in strict mode.",
    "",
    "## Rules Tested",
    "",
    "- >=1 abnormal chromosome",
    "- >=2 abnormal chromosomes",
    "- >=3 abnormal chromosomes",
    "- mean absolute log2like > threshold",
    "- sum absolute log2like > threshold",
    "",
    "## Baseline",
    "",
    sprintf("- Baseline rule (`>=1 abnormal chromosome`) sensitivity: %.3f", baseline$sensitivity[[1]]),
    sprintf("- Baseline rule (`>=1 abnormal chromosome`) specificity: %.3f", baseline$specificity[[1]]),
    sprintf("- Baseline rule (`>=1 abnormal chromosome`) accuracy: %.3f", baseline$accuracy[[1]]),
    "",
    "## Candidate Improvements",
    ""
  )

  if (nrow(improved_rules) == 0) {
    report_lines <- c(
      report_lines,
      "- No tested rule improved specificity without a substantial sensitivity drop relative to the baseline."
    )
  } else {
    for (i in seq_len(nrow(improved_rules))) {
      report_lines <- c(
        report_lines,
        sprintf(
          "- `%s`: sensitivity %.3f, specificity %.3f, accuracy %.3f",
          improved_rules$rule_name[[i]],
          improved_rules$sensitivity[[i]],
          improved_rules$specificity[[i]],
          improved_rules$accuracy[[i]]
        )
      )
    }
  }

  report_lines <- c(
    report_lines,
    "",
    "## Interpretation",
    "",
    "- Rules with higher chromosome-count requirements usually trade sensitivity for specificity.",
    "- Aggregate signal rules based on embryo-wide mean or sum absolute log2like can be more conservative than the `>=1 chromosome` rule.",
    "- Use `decision_rule_comparison.csv` for the ranked metrics and `decision_rule_confusion_matrices.csv` for the underlying confusion tables."
  )

  writeLines(report_lines, file.path(results_dir, "decision_rule_report.md"))

  message("WGSPGT decision-rule comparison complete.")
  message("  Metrics: ", file.path(results_dir, "decision_rule_comparison.csv"))
  message("  Report: ", file.path(results_dir, "decision_rule_report.md"))

  invisible(
    list(
      metrics_all = metrics_all,
      confusion_all = confusion_all
    )
  )
}

if (sys.nframe() == 0L) {
  run_wgspgt_decision_rules()
}
