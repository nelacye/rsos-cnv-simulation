# Transform workbook-derived chromosome values into a pseudo log2-like scale.
#
# Transformation:
#   1. Compute the median autosomal value (chr1-chr22) per embryo
#   2. Divide each chromosome value by that embryo-specific autosomal median
#   3. Take log2 of the normalized value

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

source("R/load_wgspgt_data.R")

run_wgspgt_log2like <- function(
    cleaned_path = file.path("results", "wgspgt_validation", "wgspgt_cleaned.csv"),
    workbook_path = file.path("data", "wgspgt", "Fig.1d_Fig.1e_Fig.1f_Fig.1g_Fig.2b_Fig5.xlsx"),
    results_dir = file.path("results", "wgspgt_descriptive")) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(cleaned_path)) {
    dat <- read_wgspgt_cleaned_data(cleaned_path)
  } else {
    dat <- save_wgspgt_cleaned_data(
      workbook_path = workbook_path,
      output_path = cleaned_path
    )
  }

  autosome_levels <- paste0("chr", 1:22)

  autosome_reference <- dat %>%
    filter(chromosome %in% autosome_levels) %>%
    group_by(embryo_id) %>%
    summarise(
      autosome_median_value = median(log2r, na.rm = TRUE),
      .groups = "drop"
    )

  transformed <- dat %>%
    left_join(autosome_reference, by = "embryo_id") %>%
    mutate(
      normalized_to_autosome_median = log2r / autosome_median_value,
      log2like = log2(normalized_to_autosome_median),
      chromosome_group = case_when(
        chromosome %in% autosome_levels ~ "autosome",
        chromosome == "chrX" ~ "chrX",
        chromosome == "chrAut" ~ "chrAut",
        TRUE ~ "other"
      )
    ) %>%
    arrange(as.numeric(embryo_id), chromosome)

  overall_summary <- transformed %>%
    summarise(
      n_rows = n(),
      min_log2like = min(log2like, na.rm = TRUE),
      q1_log2like = quantile(log2like, 0.25, na.rm = TRUE),
      median_log2like = median(log2like, na.rm = TRUE),
      mean_log2like = mean(log2like, na.rm = TRUE),
      q3_log2like = quantile(log2like, 0.75, na.rm = TRUE),
      max_log2like = max(log2like, na.rm = TRUE),
      sd_log2like = sd(log2like, na.rm = TRUE)
    )

  by_chromosome <- transformed %>%
    group_by(chromosome, chromosome_group) %>%
    summarise(
      n = n(),
      median_log2like = median(log2like, na.rm = TRUE),
      mean_log2like = mean(log2like, na.rm = TRUE),
      q1_log2like = quantile(log2like, 0.25, na.rm = TRUE),
      q3_log2like = quantile(log2like, 0.75, na.rm = TRUE),
      min_log2like = min(log2like, na.rm = TRUE),
      max_log2like = max(log2like, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(match(chromosome, c(autosome_levels, "chrX", "chrAut")))

  group_summary <- transformed %>%
    filter(chromosome_group %in% c("autosome", "chrX", "chrAut")) %>%
    group_by(chromosome_group) %>%
    summarise(
      median_log2like = median(log2like, na.rm = TRUE),
      mean_log2like = mean(log2like, na.rm = TRUE),
      q1_log2like = quantile(log2like, 0.25, na.rm = TRUE),
      q3_log2like = quantile(log2like, 0.75, na.rm = TRUE),
      min_log2like = min(log2like, na.rm = TRUE),
      max_log2like = max(log2like, na.rm = TRUE),
      .groups = "drop"
    )

  autosome_center <- group_summary %>% filter(chromosome_group == "autosome")
  chrx_center <- group_summary %>% filter(chromosome_group == "chrX")
  chraut_center <- group_summary %>% filter(chromosome_group == "chrAut")

  overall_hist <- ggplot(transformed, aes(x = log2like)) +
    geom_histogram(bins = 50, fill = "#5B7C99", color = "white") +
    theme_minimal(base_size = 11) +
    labs(
      title = "Pseudo log2-like transformed chromosome values",
      subtitle = "Overall histogram",
      x = "log2(value / embryo autosomal median)",
      y = "Count"
    )

  overall_density <- ggplot(transformed, aes(x = log2like, color = chromosome_group, fill = chromosome_group)) +
    geom_density(alpha = 0.2) +
    theme_minimal(base_size = 11) +
    scale_color_manual(values = c("autosome" = "#406882", "chrX" = "#9A031E", "chrAut" = "#5F0F40", "other" = "#999999")) +
    scale_fill_manual(values = c("autosome" = "#406882", "chrX" = "#9A031E", "chrAut" = "#5F0F40", "other" = "#999999")) +
    labs(
      title = "Pseudo log2-like transformed values by chromosome group",
      x = "log2(value / embryo autosomal median)",
      y = "Density",
      color = "Group",
      fill = "Group"
    )

  chromosome_boxplot <- ggplot(transformed, aes(x = reorder(chromosome, log2like, median), y = log2like, fill = chromosome_group)) +
    geom_boxplot(outlier.alpha = 0.35) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = c("autosome" = "#7AA6C2", "chrX" = "#D1495B", "chrAut" = "#6D597A", "other" = "#999999")) +
    labs(
      title = "Pseudo log2-like values by chromosome",
      x = "Chromosome",
      y = "log2(value / embryo autosomal median)",
      fill = "Group"
    )

  ggsave(file.path(results_dir, "wgspgt_log2like_histogram.png"), overall_hist, width = 8, height = 5, dpi = 300)
  ggsave(file.path(results_dir, "wgspgt_log2like_density.png"), overall_density, width = 8, height = 5, dpi = 300)
  ggsave(file.path(results_dir, "wgspgt_log2like_by_chromosome.png"), chromosome_boxplot, width = 11, height = 5, dpi = 300)

  output_transformed <- transformed %>%
    select(
      embryo_id,
      chromosome,
      log2r,
      autosome_median_value,
      normalized_to_autosome_median,
      log2like,
      meanDepth,
      meanBreadth,
      EmbryoStatus,
      Method
    )

  write_csv(output_transformed, file.path(results_dir, "wgspgt_log2like.csv"))
  write_csv(by_chromosome, file.path(results_dir, "wgspgt_log2like_by_chromosome.csv"))
  write_csv(group_summary, file.path(results_dir, "wgspgt_log2like_group_summary.csv"))

  autosome_median_abs <- median(abs(transformed$log2like[transformed$chromosome_group == "autosome"]), na.rm = TRUE)
  chrx_median <- chrx_center$median_log2like[[1]]
  chraut_median <- chraut_center$median_log2like[[1]]

  plausibility_lines <- c(
    "# Pseudo log2-like Transformation Report",
    "",
    "## Transformation",
    "",
    "For each embryo:",
    "",
    "1. Compute the median autosomal value across `chr1`-`chr22`.",
    "2. Divide each chromosome value by that embryo-specific autosomal median.",
    "3. Take `log2` of the normalized value.",
    "",
    "## Distribution Check",
    "",
    sprintf("- Overall median transformed value: %.3f", overall_summary$median_log2like[[1]]),
    sprintf("- Overall interquartile range: %.3f to %.3f", overall_summary$q1_log2like[[1]], overall_summary$q3_log2like[[1]]),
    sprintf("- Autosome group median: %.3f", autosome_center$median_log2like[[1]]),
    sprintf("- Autosome median absolute deviation from zero (using absolute transformed values): %.3f", autosome_median_abs),
    sprintf("- `chrAut` median transformed value: %.3f", chraut_median),
    sprintf("- `chrX` median transformed value: %.3f", chrx_median),
    "",
    "## Interpretation",
    "",
    if (abs(autosome_center$median_log2like[[1]]) < 0.05) {
      "- Autosomes are centered very close to 0 after transformation, which is consistent with a log2-ratio-like baseline."
    } else {
      "- Autosomes are not centered close enough to 0 to cleanly mimic a standard log2-ratio baseline."
    },
    if (abs(chraut_median) < 0.15) {
      "- `chrAut` remains close to the autosomal baseline after transformation."
    } else {
      "- `chrAut` remains shifted relative to the autosomal baseline after transformation."
    },
    if (abs(chrx_median) < 1.5) {
      "- `chrX` moves into a more interpretable range after transformation."
    } else {
      "- `chrX` remains far from a conventional CNV-like range after transformation."
    },
    "- This normalization improves comparability across embryos because it removes embryo-specific scale offsets.",
    "- However, the transformed values should still be treated as pseudo log2-like rather than validated QDNAseq logR values.",
    "",
    "## Practical Conclusion",
    "",
    if (abs(autosome_center$median_log2like[[1]]) < 0.05 && abs(chrx_median) < 1.5) {
      "- The transformation yields a plausible first-pass log2-like signal for exploratory downstream CNV work."
    } else {
      "- The transformation partially improves the scale, but it does not fully establish a trustworthy log2-like CNV signal."
    },
    "- Use `wgspgt_log2like.csv` for exploratory modeling only until the original QDNAseq `LogR` or `SegLogR` outputs can be recovered."
  )

  writeLines(plausibility_lines, file.path(results_dir, "wgspgt_log2like_report.md"))

  message("WGS-PGT pseudo log2-like transformation complete.")
  message("  Output dataset: ", file.path(results_dir, "wgspgt_log2like.csv"))
  message("  Report: ", file.path(results_dir, "wgspgt_log2like_report.md"))

  invisible(
    list(
      transformed = output_transformed,
      by_chromosome = by_chromosome,
      group_summary = group_summary
    )
  )
}

if (sys.nframe() == 0L) {
  run_wgspgt_log2like()
}
