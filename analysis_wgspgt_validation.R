# First-pass validation workflow for the WGS-PGT workbook.
#
# This analysis intentionally skips diagnosis-sheet logic for now. It creates a
# cleaned chromosome-level dataset, applies the existing threshold classifier,
# and saves descriptive outputs that prepare a later Bayesian-vs-threshold
# comparison.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("R/classify_embryo.R")
source("R/load_wgspgt_data.R")

run_wgspgt_validation <- function(
    workbook_path = file.path("data", "wgspgt", "Fig.1d_Fig.1e_Fig.1f_Fig.1g_Fig.2b_Fig5.xlsx"),
    low_thr = 0.15,
    high_thr = 0.30,
    results_dir = file.path("results", "wgspgt_validation")) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  cleaned_path <- file.path(results_dir, "wgspgt_cleaned.csv")

  # Build the cleaned validation dataset from the two relevant workbook sheets.
  cleaned_data <- save_wgspgt_cleaned_data(
    workbook_path = workbook_path,
    output_path = cleaned_path
  )

  # Read back the cleaned data so downstream steps run on the stable CSV.
  validation_data <- read_wgspgt_cleaned_data(cleaned_path)

  threshold_calls <- validation_data %>%
    mutate(
      threshold_call = classify_embryo(log2r, low_thr = low_thr, high_thr = high_thr)
    )

  chromosome_summary <- threshold_calls %>%
    count(Method, EmbryoStatus, threshold_call, name = "n_chromosomes") %>%
    arrange(Method, EmbryoStatus, threshold_call)

  embryo_summary <- threshold_calls %>%
    group_by(embryo_id, Method, EmbryoStatus, meanDepth, meanBreadth) %>%
    summarise(
      n_chromosomes = n(),
      mean_abs_log2r = mean(abs_log2r, na.rm = TRUE),
      max_abs_log2r = max(abs_log2r, na.rm = TRUE),
      n_normal = sum(threshold_call == "normal", na.rm = TRUE),
      n_mosaic = sum(threshold_call == "mosaic", na.rm = TRUE),
      n_aneuploid = sum(threshold_call == "aneuploid", na.rm = TRUE),
      any_abnormal_threshold = any(threshold_call != "normal", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(embryo_id)

  method_summary <- embryo_summary %>%
    group_by(Method, EmbryoStatus) %>%
    summarise(
      n_embryos = n(),
      mean_depth = mean(meanDepth, na.rm = TRUE),
      mean_breadth = mean(meanBreadth, na.rm = TRUE),
      mean_abs_log2r = mean(mean_abs_log2r, na.rm = TRUE),
      embryos_with_any_abnormal_threshold = sum(any_abnormal_threshold, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Method, EmbryoStatus)

  readr::write_csv(threshold_calls, file.path(results_dir, "wgspgt_threshold_calls.csv"))
  readr::write_csv(chromosome_summary, file.path(results_dir, "wgspgt_chromosome_summary.csv"))
  readr::write_csv(embryo_summary, file.path(results_dir, "wgspgt_embryo_summary.csv"))
  readr::write_csv(method_summary, file.path(results_dir, "wgspgt_method_summary.csv"))

  message("WGS-PGT first-pass validation complete.")
  message("  Cleaned data: ", cleaned_path)
  message("  Threshold calls: ", file.path(results_dir, "wgspgt_threshold_calls.csv"))
  message("  Chromosome summary: ", file.path(results_dir, "wgspgt_chromosome_summary.csv"))
  message("  Embryo summary: ", file.path(results_dir, "wgspgt_embryo_summary.csv"))
  message("  Method summary: ", file.path(results_dir, "wgspgt_method_summary.csv"))

  invisible(
    list(
      cleaned_data = validation_data,
      threshold_calls = threshold_calls,
      chromosome_summary = chromosome_summary,
      embryo_summary = embryo_summary,
      method_summary = method_summary
    )
  )
}

args <- commandArgs(trailingOnly = TRUE)

if (sys.nframe() == 0L) {
  workbook_arg <- if (length(args) >= 1) {
    args[[1]]
  } else {
    file.path("data", "wgspgt", "Fig.1d_Fig.1e_Fig.1f_Fig.1g_Fig.2b_Fig5.xlsx")
  }

  run_wgspgt_validation(workbook_path = workbook_arg)
}
