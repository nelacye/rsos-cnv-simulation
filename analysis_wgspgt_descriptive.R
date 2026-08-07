# Descriptive analysis of the workbook-specific WGS-PGT chromosome values.
#
# This script uses the cleaned dataset derived from:
#   - Fig1f_WGS
#   - Fig1d_Fig1e
# and intentionally does not apply any classifier yet.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

source("R/load_wgspgt_data.R")

safe_density_path <- function(results_dir, filename) {
  file.path(results_dir, filename)
}

run_wgspgt_descriptive <- function(
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

  dat <- dat %>%
    mutate(
      chromosome_group = case_when(
        chromosome %in% autosome_levels ~ "autosome",
        chromosome == "chrX" ~ "chrX",
        chromosome == "chrAut" ~ "chrAut",
        TRUE ~ "other"
      ),
      EmbryoStatus = if_else(is.na(EmbryoStatus), "missing", EmbryoStatus),
      Method = if_else(is.na(Method), "missing", Method)
    )

  overall_summary <- dat %>%
    summarise(
      n_rows = n(),
      n_embryos = n_distinct(embryo_id),
      min_value = min(log2r, na.rm = TRUE),
      q1_value = quantile(log2r, 0.25, na.rm = TRUE),
      median_value = median(log2r, na.rm = TRUE),
      mean_value = mean(log2r, na.rm = TRUE),
      q3_value = quantile(log2r, 0.75, na.rm = TRUE),
      max_value = max(log2r, na.rm = TRUE),
      sd_value = sd(log2r, na.rm = TRUE)
    )

  by_chromosome <- dat %>%
    group_by(chromosome, chromosome_group) %>%
    summarise(
      n = n(),
      min_value = min(log2r, na.rm = TRUE),
      q1_value = quantile(log2r, 0.25, na.rm = TRUE),
      median_value = median(log2r, na.rm = TRUE),
      mean_value = mean(log2r, na.rm = TRUE),
      q3_value = quantile(log2r, 0.75, na.rm = TRUE),
      max_value = max(log2r, na.rm = TRUE),
      sd_value = sd(log2r, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(match(chromosome, c(autosome_levels, "chrX", "chrAut")))

  by_status_method <- dat %>%
    group_by(EmbryoStatus, Method) %>%
    summarise(
      n = n(),
      n_embryos = n_distinct(embryo_id),
      median_value = median(log2r, na.rm = TRUE),
      mean_value = mean(log2r, na.rm = TRUE),
      q1_value = quantile(log2r, 0.25, na.rm = TRUE),
      q3_value = quantile(log2r, 0.75, na.rm = TRUE),
      max_value = max(log2r, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Method, EmbryoStatus)

  by_status_chromosome_group <- dat %>%
    group_by(EmbryoStatus, chromosome_group) %>%
    summarise(
      n = n(),
      median_value = median(log2r, na.rm = TRUE),
      mean_value = mean(log2r, na.rm = TRUE),
      q1_value = quantile(log2r, 0.25, na.rm = TRUE),
      q3_value = quantile(log2r, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(EmbryoStatus, chromosome_group)

  embryo_level_summary <- dat %>%
    group_by(embryo_id, EmbryoStatus, Method, meanDepth, meanBreadth) %>%
    mutate(
      embryo_q1 = quantile(log2r, 0.25, na.rm = TRUE),
      embryo_q3 = quantile(log2r, 0.75, na.rm = TRUE),
      embryo_iqr = embryo_q3 - embryo_q1,
      low_fence = embryo_q1 - 1.5 * embryo_iqr,
      high_fence = embryo_q3 + 1.5 * embryo_iqr,
      is_extreme_chromosome = log2r < low_fence | log2r > high_fence
    ) %>%
    summarise(
      n_chromosomes = n(),
      median_chromosome_value = median(log2r, na.rm = TRUE),
      min_chromosome_value = min(log2r, na.rm = TRUE),
      max_chromosome_value = max(log2r, na.rm = TRUE),
      chromosome_value_range = max(log2r, na.rm = TRUE) - min(log2r, na.rm = TRUE),
      n_extreme_chromosomes = sum(is_extreme_chromosome, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(as.numeric(embryo_id))

  autosome_reference <- dat %>%
    filter(chromosome_group == "autosome") %>%
    summarise(
      autosome_median = median(log2r, na.rm = TRUE),
      autosome_q1 = quantile(log2r, 0.25, na.rm = TRUE),
      autosome_q3 = quantile(log2r, 0.75, na.rm = TRUE)
    )

  scale_comparison <- dat %>%
    filter(chromosome_group %in% c("autosome", "chrX", "chrAut")) %>%
    group_by(chromosome_group) %>%
    summarise(
      median_value = median(log2r, na.rm = TRUE),
      mean_value = mean(log2r, na.rm = TRUE),
      q1_value = quantile(log2r, 0.25, na.rm = TRUE),
      q3_value = quantile(log2r, 0.75, na.rm = TRUE),
      max_value = max(log2r, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      autosome_median = autosome_reference$autosome_median,
      median_ratio_to_autosome = median_value / autosome_reference$autosome_median
    )

  top_extreme_embryos <- embryo_level_summary %>%
    arrange(desc(n_extreme_chromosomes), desc(chromosome_value_range), embryo_id)

  overall_hist <- ggplot(dat, aes(x = log2r)) +
    geom_histogram(bins = 50, fill = "#5B7C99", color = "white") +
    theme_minimal(base_size = 11) +
    labs(
      title = "Fig1f_WGS chromosome-level values",
      subtitle = "Overall histogram",
      x = "Chromosome-level value",
      y = "Count"
    )

  overall_density <- ggplot(dat, aes(x = log2r)) +
    geom_density(fill = "#97B8D1", alpha = 0.7, color = "#335C81") +
    theme_minimal(base_size = 11) +
    labs(
      title = "Fig1f_WGS chromosome-level values",
      subtitle = "Overall density",
      x = "Chromosome-level value",
      y = "Density"
    )

  chromosome_density <- ggplot(dat, aes(x = log2r, group = chromosome_group, color = chromosome_group, fill = chromosome_group)) +
    geom_density(alpha = 0.2) +
    theme_minimal(base_size = 11) +
    scale_color_manual(values = c("autosome" = "#406882", "chrX" = "#9A031E", "chrAut" = "#5F0F40", "other" = "#999999")) +
    scale_fill_manual(values = c("autosome" = "#406882", "chrX" = "#9A031E", "chrAut" = "#5F0F40", "other" = "#999999")) +
    labs(
      title = "Chromosome-group density comparison",
      x = "Chromosome-level value",
      y = "Density",
      color = "Group",
      fill = "Group"
    )

  by_chromosome_boxplot <- ggplot(dat, aes(x = reorder(chromosome, log2r, median), y = log2r, fill = chromosome_group)) +
    geom_boxplot(outlier.alpha = 0.4) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = c("autosome" = "#7AA6C2", "chrX" = "#D1495B", "chrAut" = "#6D597A", "other" = "#999999")) +
    labs(
      title = "Chromosome-level values by chromosome",
      x = "Chromosome",
      y = "Chromosome-level value",
      fill = "Group"
    )

  by_status_density <- ggplot(dat, aes(x = log2r, color = EmbryoStatus, fill = EmbryoStatus)) +
    geom_density(alpha = 0.2) +
    theme_minimal(base_size = 11) +
    labs(
      title = "Chromosome-level values by EmbryoStatus",
      x = "Chromosome-level value",
      y = "Density"
    )

  by_status_boxplot <- ggplot(dat, aes(x = EmbryoStatus, y = log2r, fill = EmbryoStatus)) +
    geom_boxplot(outlier.alpha = 0.4) +
    theme_minimal(base_size = 11) +
    labs(
      title = "Chromosome-level values by EmbryoStatus",
      x = "EmbryoStatus",
      y = "Chromosome-level value"
    )

  embryo_median_plot <- ggplot(embryo_level_summary, aes(x = reorder(embryo_id, median_chromosome_value), y = median_chromosome_value, fill = EmbryoStatus)) +
    geom_col() +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
      title = "Embryo-level median chromosome value",
      x = "Embryo ID",
      y = "Median chromosome value"
    )

  embryo_range_plot <- ggplot(embryo_level_summary, aes(x = reorder(embryo_id, chromosome_value_range), y = chromosome_value_range, fill = EmbryoStatus)) +
    geom_col() +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
      title = "Embryo-level chromosome value range",
      x = "Embryo ID",
      y = "Range across chromosomes"
    )

  embryo_extreme_plot <- ggplot(embryo_level_summary, aes(x = reorder(embryo_id, n_extreme_chromosomes), y = n_extreme_chromosomes, fill = EmbryoStatus)) +
    geom_col() +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
      title = "Embryo-level extreme chromosome counts",
      subtitle = "Extreme defined by embryo-specific Tukey fences",
      x = "Embryo ID",
      y = "Number of extreme chromosomes"
    )

  ggsave(safe_density_path(results_dir, "overall_histogram.png"), overall_hist, width = 8, height = 5, dpi = 300)
  ggsave(safe_density_path(results_dir, "overall_density.png"), overall_density, width = 8, height = 5, dpi = 300)
  ggsave(safe_density_path(results_dir, "chromosome_group_density.png"), chromosome_density, width = 8, height = 5, dpi = 300)
  ggsave(safe_density_path(results_dir, "chromosome_boxplot.png"), by_chromosome_boxplot, width = 11, height = 5, dpi = 300)
  ggsave(safe_density_path(results_dir, "status_density.png"), by_status_density, width = 8, height = 5, dpi = 300)
  ggsave(safe_density_path(results_dir, "status_boxplot.png"), by_status_boxplot, width = 7, height = 5, dpi = 300)
  ggsave(safe_density_path(results_dir, "embryo_median_values.png"), embryo_median_plot, width = 10, height = 5, dpi = 300)
  ggsave(safe_density_path(results_dir, "embryo_value_ranges.png"), embryo_range_plot, width = 10, height = 5, dpi = 300)
  ggsave(safe_density_path(results_dir, "embryo_extreme_counts.png"), embryo_extreme_plot, width = 10, height = 5, dpi = 300)

  readr::write_csv(overall_summary, file.path(results_dir, "overall_summary.csv"))
  readr::write_csv(by_chromosome, file.path(results_dir, "by_chromosome_summary.csv"))
  readr::write_csv(by_status_method, file.path(results_dir, "by_status_method_summary.csv"))
  readr::write_csv(by_status_chromosome_group, file.path(results_dir, "by_status_chromosome_group_summary.csv"))
  readr::write_csv(embryo_level_summary, file.path(results_dir, "embryo_level_summary.csv"))
  readr::write_csv(scale_comparison, file.path(results_dir, "chromosome_scale_comparison.csv"))

  method_count <- n_distinct(dat$Method)
  embryo_status_levels <- paste(sort(unique(dat$EmbryoStatus)), collapse = ", ")
  chrx_ratio <- scale_comparison %>% filter(chromosome_group == "chrX") %>% pull(median_ratio_to_autosome)
  chraut_ratio <- scale_comparison %>% filter(chromosome_group == "chrAut") %>% pull(median_ratio_to_autosome)
  overall_median <- overall_summary$median_value[[1]]
  overall_max <- overall_summary$max_value[[1]]
  autosome_median_min <- min(by_chromosome$median_value[by_chromosome$chromosome_group == "autosome"])
  autosome_median_max <- max(by_chromosome$median_value[by_chromosome$chromosome_group == "autosome"])
  top_embryo <- top_extreme_embryos[1, ]

  report_lines <- c(
    "# WGS-PGT Descriptive Report",
    "",
    "## Scope",
    "",
    "This report uses the workbook-specific cleaned dataset derived from `Fig1f_WGS` and `Fig1d_Fig1e` only. No Bayesian classifier or threshold classifier was applied in this analysis.",
    "",
    "Outputs in this directory:",
    "",
    "- `overall_summary.csv`",
    "- `by_chromosome_summary.csv`",
    "- `by_status_method_summary.csv`",
    "- `by_status_chromosome_group_summary.csv`",
    "- `embryo_level_summary.csv`",
    "- `chromosome_scale_comparison.csv`",
    "- plot PNG files",
    "",
    "## Dataset Overview",
    "",
    sprintf("- Chromosome-level rows: %s", comma(overall_summary$n_rows[[1]])),
    sprintf("- Embryos: %s", comma(overall_summary$n_embryos[[1]])),
    sprintf("- EmbryoStatus levels: %s", embryo_status_levels),
    sprintf("- Method levels observed: %s", method_count),
    "",
    "## Overall Distribution",
    "",
    sprintf("- Median value: %.3f", overall_median),
    sprintf("- Mean value: %.3f", overall_summary$mean_value[[1]]),
    sprintf("- Interquartile range: %.3f to %.3f", overall_summary$q1_value[[1]], overall_summary$q3_value[[1]]),
    sprintf("- Minimum to maximum: %.3f to %.3f", overall_summary$min_value[[1]], overall_max),
    "- The distribution is strongly right-skewed with a long upper tail.",
    "",
    "## By Chromosome",
    "",
    sprintf("- Autosomes are broadly concentrated around medians of roughly %.1f to %.1f.", autosome_median_min, autosome_median_max),
    sprintf("- `chrX` median is %.3f, which is %.1fx the pooled autosome median.", scale_comparison$median_value[scale_comparison$chromosome_group == "chrX"], chrx_ratio),
    sprintf("- `chrAut` median is %.3f, which is %.1fx the pooled autosome median.", scale_comparison$median_value[scale_comparison$chromosome_group == "chrAut"], chraut_ratio),
    "- `chrX` is clearly on a different scale from autosomes.",
    "- `chrAut` is elevated relative to autosomes but remains much closer to the autosome range than `chrX`.",
    "",
    "## By EmbryoStatus and Method",
    "",
    if (method_count == 1) {
      "- Only one method (`WGS`) is present, so there is no between-method distribution comparison available in this dataset."
    } else {
      "- Multiple methods are present; see `by_status_method_summary.csv` for method-level comparisons."
    },
    "- Distribution differences across `EmbryoStatus` are visible, but they occur on top of the strong chromosome-specific scale effects, especially `chrX`.",
    "",
    "## Embryo-Level Summaries",
    "",
    "- For each embryo, the report includes:",
    "  - median chromosome value",
    "  - range across chromosomes",
    "  - number of extreme chromosomes",
    "- Extreme chromosomes were defined using embryo-specific Tukey fences: values below Q1 - 1.5*IQR or above Q3 + 1.5*IQR within that embryo.",
    sprintf("- The embryo with the most extreme chromosomes in this first pass is embryo `%s` with %s extreme chromosomes and a chromosome-value range of %.3f.", top_embryo$embryo_id[[1]], top_embryo$n_extreme_chromosomes[[1]], top_embryo$chromosome_value_range[[1]]),
    "",
    "## Interpretation of Fig1f_WGS Scale",
    "",
    "- These values do not behave like standard log2 ratios.",
    sprintf("  Standard CNV log2 ratios are typically centered near 0 for diploid regions, whereas these values have a pooled median of %.3f and extend to %.3f.", overall_median, overall_max),
    "- The autosomes are also not centered near a simple diploid baseline such as 0 or 2-copy-equivalent log2 space.",
    "- The `chrX` values are dramatically inflated relative to autosomes, which argues against this being a uniformly normalized log2-ratio track.",
    "- The data do show chromosome-structured behavior rather than looking like a pure random QC/error score, so they are unlikely to be a generic sequencing QC metric alone.",
    "- The most plausible interpretation from this workbook alone is that `Fig1f_WGS` behaves more like a copy-number-like or coverage-derived signal on a transformed positive scale, not a standard log2-ratio output.",
    "- Because the scale is positive, shifted, and chromosome-dependent, it should not be fed directly into a Bayesian CNV classifier that expects conventional log2 ratios without a prior scale-reconstruction step.",
    "",
    "## Practical Conclusion",
    "",
    "- `Fig1f_WGS` is useful for descriptive exploration and embryo/chromosome ranking.",
    "- It is not yet suitable as a drop-in replacement for a standard chromosome-level log2-ratio input.",
    "- Before Bayesian-vs-threshold comparison, the next technical step should be to recover or reconstruct the original QDNAseq `LogR`/`SegLogR` scale, or at minimum derive a validated mapping from this workbook scale to a CNV-relevant measurement."
  )

  writeLines(report_lines, file.path(results_dir, "wgspgt_descriptive_report.md"))

  message("WGS-PGT descriptive analysis complete.")
  message("  Results directory: ", results_dir)
  message("  Report: ", file.path(results_dir, "wgspgt_descriptive_report.md"))

  invisible(
    list(
      overall_summary = overall_summary,
      by_chromosome = by_chromosome,
      by_status_method = by_status_method,
      by_status_chromosome_group = by_status_chromosome_group,
      embryo_level_summary = embryo_level_summary,
      scale_comparison = scale_comparison
    )
  )
}

if (sys.nframe() == 0L) {
  run_wgspgt_descriptive()
}
