# ============================================================
# Borderline CNVs in PGT-A: A Simulation Study
# Royal Society Open Science
# 
# This script simulates the effect of current threshold-based
# classification on mosaic CNV detection in preimplantation
# genetic testing.
# 
# Author: Anere Cye
# Reproducible at: https://github.com/anerecye/rsos-cnv-simulation
# ============================================================

# 1. Setup -----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)

# Create output directories
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

set.seed(2025)

# 2. Simulation function ---------------------------------------------------

#' Simulate observed log2 ratio for a single embryo
#'
#' @param true_mosaicism Fraction of cells carrying the CNV (0 to 1)
#' @param cnv_type Either "dup" (duplication) or "del" (deletion)
#' @param coverage Mean sequencing coverage (e.g., 10x)
#' @param noise_sd_base Baseline standard deviation at 10x coverage
#' @param pcr_duplicates Multiplier for PCR-induced noise
#' @return Observed log2 ratio after adding technical noise

simulate_cnv <- function(true_mosaicism, cnv_type = "dup",
                         coverage = 10, noise_sd_base = 0.05,
                         pcr_duplicates = 1.2) {
  
  # True copy ratio (CR) relative to diploid control
  # For duplication (+1 copy): CR = (M*3 + (1-M)*2)/2 = 1 + M/2
  # For deletion (-1 copy):  CR = (M*1 + (1-M)*2)/2 = 1 - M/2
  
  if (cnv_type == "dup") {
    true_cr <- 1 + true_mosaicism / 2
  } else if (cnv_type == "del") {
    true_cr <- 1 - true_mosaicism / 2
  } else {
    stop("cnv_type must be 'dup' or 'del'")
  }
  
  # Noise inversely proportional to sqrt(coverage)
  effective_noise_sd <- noise_sd_base * pcr_duplicates / sqrt(coverage / 10)
  
  # Add Gaussian noise
  observed_cr <- true_cr + rnorm(1, mean = 0, sd = effective_noise_sd)
  
  # Log2 transformation for threshold application
  log2_ratio <- log2(observed_cr)
  
  return(log2_ratio)
}

# 3. Experimental design ---------------------------------------------------

# True mosaicism levels (0% to 30%)
mosaicism_levels <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

# Sequencing coverages (low to high, realistic for PGT-A)
coverages <- c(5, 10, 20, 50)

# CNV type (focus on duplications; deletions symmetric)
cnv_types <- c("dup")

# Classification thresholds (lower and upper bounds)
low_thresholds <- c(0.10, 0.15, 0.20)
high_thresholds <- c(0.30, 0.40)

# Number of replicates per parameter combination
n_replicates <- 5000

# 4. Run simulation --------------------------------------------------------

# Create all combinations
simulation_grid <- expand_grid(
  mosaicism = mosaicism_levels,
  coverage = coverages,
  rep = 1:n_replicates,
  cnv_type = cnv_types
)

# Apply simulation to each row
sim_results <- simulation_grid %>%
  rowwise() %>%
  mutate(
    log2r = simulate_cnv(
      true_mosaicism = mosaicism,
      cnv_type = cnv_type,
      coverage = coverage
    )
  ) %>%
  ungroup()

# Save raw results (optional, can be large ~2GB)
write_csv(sim_results, "results/simulated_log2ratios.csv")

# 5. Apply classification thresholds ---------------------------------------

# Create all threshold combinations
threshold_pairs <- expand_grid(
  low = low_thresholds,
  high = high_thresholds
) %>%
  filter(low < high)

# Classify function
classify_embryo <- function(log2r, low_thr, high_thr) {
  case_when(
    abs(log2r) < low_thr   ~ "normal",
    abs(log2r) < high_thr  ~ "mosaic",
    TRUE                   ~ "aneuploid"
  )
}

# Apply to all threshold pairs
all_classifications <- bind_rows(lapply(1:nrow(threshold_pairs), function(i) {
  low_thr <- threshold_pairs$low[i]
  high_thr <- threshold_pairs$high[i]
  
  sim_results %>%
    mutate(
      class = classify_embryo(log2r, low_thr, high_thr),
      low_threshold = low_thr,
      high_threshold = high_thr
    )
}))

# 6. Calculate metrics with confidence intervals --------------------------

# Function for Wilson confidence interval
wilson_ci <- function(n_success, n_total, conf_level = 0.95) {
  prop.test(n_success, n_total, conf.level = conf_level)$conf.int
}

# Aggregate proportions for each class
compute_proportions <- function(df, class_label) {
  df %>%
    group_by(coverage, low_threshold, high_threshold, mosaicism) %>%
    summarise(
      n_total = n(),
      n_event = sum(class == class_label),
      proportion = n_event / n_total,
      ci_lower = wilson_ci(n_event, n_total)[1],
      ci_upper = wilson_ci(n_event, n_total)[2],
      .groups = "drop"
    ) %>%
    mutate(class = class_label)
}

prop_normal <- compute_proportions(all_classifications, "normal")
prop_mosaic <- compute_proportions(all_classifications, "mosaic")
prop_aneuploid <- compute_proportions(all_classifications, "aneuploid")

# Combine into wide format
full_metrics <- bind_rows(prop_normal, prop_mosaic, prop_aneuploid) %>%
  pivot_wider(
    id_cols = c(coverage, low_threshold, high_threshold, mosaicism),
    names_from = class,
    values_from = c(proportion, ci_lower, ci_upper),
    names_sep = "_"
  )

# Save full metrics
write_csv(full_metrics, "results/full_classification_metrics.csv")

# 7. Robustness analysis ---------------------------------------------------

# Define critical threshold: more than 50% misclassified as normal
critical_threshold <- 0.50

robustness_check <- full_metrics %>%
  filter(mosaicism %in% c(0.15, 0.20)) %>%
  mutate(
    is_critical = proportion_normal > critical_threshold,
    ci_above_critical = ci_lower_normal > critical_threshold
  ) %>%
  select(coverage, low_threshold, high_threshold, mosaicism,
         proportion_normal, ci_lower_normal, ci_upper_normal,
         is_critical, ci_above_critical)

# Summary by mosaicism level
robustness_summary <- robustness_check %>%
  group_by(mosaicism) %>%
  summarise(
    n_combinations = n(),
    n_ci_above_critical = sum(ci_above_critical, na.rm = TRUE),
    percent_robust = 100 * n_ci_above_critical / n_combinations,
    .groups = "drop"
  )

# Save robustness results
write_csv(robustness_check, "results/robustness_check.csv")
write_csv(robustness_summary, "results/robustness_summary.csv")

# 8. Figures for publication -----------------------------------------------

# Figure 1: Stacked bar plot for default thresholds (low=0.15, high=0.30)
default_low <- 0.15
default_high <- 0.30

fig1_data <- full_metrics %>%
  filter(low_threshold == default_low, high_threshold == default_high)

p1 <- ggplot(fig1_data, aes(x = factor(mosaicism * 100),
                            y = proportion_normal, fill = "Normal")) +
  geom_col(aes(y = proportion_normal, fill = "Normal"), width = 0.7) +
  geom_col(aes(y = proportion_mosaic, fill = "Mosaic"), width = 0.7,
           position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = proportion_aneuploid, fill = "Aneuploid"), width = 0.7,
           position = position_stack(reverse = TRUE)) +
  facet_wrap(~coverage, labeller = label_both) +
  scale_fill_manual(values = c("Normal" = "#2E8B57",
                               "Mosaic" = "#FFA500",
                               "Aneuploid" = "#CD5C5C")) +
  labs(x = "True mosaicism level (%)",
       y = "Proportion of classified embryos",
       fill = "Class",
       title = paste0("Classification outcomes at thresholds: normal < ", default_low,
                      ", mosaic < ", default_high)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray90"))

ggsave("figures/fig1_classification_stacked_bars.png", p1, width = 10, height = 6, dpi = 300)

# Figure 2: Trade-off between sensitivity and false positive rate
trade_off_data <- all_classifications %>%
  filter(high_threshold == 0.30, coverage %in% c(5, 10, 20)) %>%
  group_by(coverage, low_threshold, mosaicism) %>%
  summarise(
    sensitivity = mean(class != "normal"),
    .groups = "drop"
  ) %>%
  left_join(
    all_classifications %>%
      filter(high_threshold == 0.30, coverage %in% c(5, 10, 20), mosaicism == 0) %>%
      group_by(coverage, low_threshold) %>%
      summarise(fpr = mean(class != "normal"), .groups = "drop"),
    by = c("coverage", "low_threshold")
  ) %>%
  filter(mosaicism %in% c(0.10, 0.15, 0.20))

p2 <- ggplot(trade_off_data, aes(x = fpr, y = sensitivity,
                                 color = factor(coverage),
                                 shape = factor(mosaicism * 100))) +
  geom_line(aes(group = interaction(coverage, mosaicism)), alpha = 0.5) +
  geom_point(size = 3) +
  scale_x_continuous(labels = percent_format(), limits = c(0, 0.5)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(x = "False positive rate (among diploid embryos, M = 0%)",
       y = "Sensitivity (proportion called mosaic or aneuploid)",
       color = "Coverage (x)",
       shape = "True mosaicism (%)",
       title = "Diagnostic trade-off: lowering the threshold increases both sensitivity and FPR") +
  theme_minimal(base_size = 12)

ggsave("figures/fig2_tradeoff_curve.png", p2, width = 8, height = 5, dpi = 300)

# Figure 3: Heatmap of false negatives
fig3_data <- full_metrics %>%
  filter(low_threshold == default_low, high_threshold == default_high) %>%
  mutate(false_normal_rate = proportion_normal)

p3 <- ggplot(fig3_data, aes(x = factor(mosaicism * 100),
                            y = factor(coverage),
                            fill = false_normal_rate)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "lightcoral", high = "darkred",
                       midpoint = 0.3, labels = percent_format(),
                       name = "False normal rate") +
  labs(x = "True mosaicism level (%)",
       y = "Coverage (x)",
       title = "Risk of missing mosaic embryos (false negatives)") +
  theme_minimal(base_size = 12) +
  coord_fixed()

ggsave("figures/fig3_heatmap_false_negatives.png", p3, width = 6, height = 5, dpi = 300)

# 9. Summary table for manuscript ------------------------------------------

summary_table <- full_metrics %>%
  filter(coverage == 10,
         low_threshold == 0.15,
         high_threshold == 0.30) %>%
  select(mosaicism, proportion_normal, proportion_mosaic, proportion_aneuploid) %>%
  mutate(across(starts_with("proportion"), ~ percent(.x, accuracy = 0.1)))

write_csv(summary_table, "results/table_classification_rates.csv")

# Print to console
cat("\n========================================\n")
cat("SIMULATION COMPLETE\n")
cat("========================================\n")
cat("Output files saved in 'figures/' and 'results/'\n\n")

cat("Robustness summary:\n")
print(robustness_summary)

cat("\n All done. Ready for Royal Society Open Science submission.\n")
