
R version 4.5.3 (2026-03-11 ucrt) -- "Reassured Reassurer"
Copyright (C) 2026 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> # ============================================================
> # Sensitivity Playground
> # 
> # Modify the parameters below to explore how they affect
> # classification outcomes. No changes to main code needed.
> # 
> # Run this script line by line or source() it.
> # ============================================================
> # ============================================
> # USER-ADJUSTABLE PARAMETERS
> # ============================================
> 
> # Noise parameters
> noise_sd_base <- 0.05      # try: 0.03, 0.04, 0.06, 0.08
> pcr_duplicates <- 1.2      # try: 1.0, 1.3, 1.5
> batch_sd <- 0.05           # try: 0, 0.02, 0.10
> correlated_noise <- TRUE   # try: FALSE
> 
> # Technical parameters
> coverage <- 10             # try: 5, 10, 20, 50
> n_replicates <- 1000       # try: 500, 2000, 5000 (higher = slower)
> 
> # Thresholds
> low_thr <- 0.15            # try: 0.10, 0.12, 0.18, 0.20
> high_thr <- 0.30           # try: 0.25, 0.35, 0.40
> 
> # Scenario selection (which one to run)
> run_single_chr <- TRUE
> run_scenario_B <- TRUE
> 
> # ============================================
> # ============================================
> # LOAD FUNCTIONS
> # ============================================
> 
> source("R/simulate_cnv.R")
> source("R/classify_embryo.R")
> source("R/simulate_embryo_single_affected.R")
> 
> set.seed(2025)
> # ============================================
> # SINGLE CHROMOSOME SIMULATION (custom noise)
> # ============================================
> 
> simulate_cnv_custom <- function(true_mosaicism, cnv_type = "dup") {
+   if (cnv_type == "dup") {
+     true_cr <- 1 + true_mosaicism / 2
+   } else {
+     true_cr <- 1 - true_mosaicism / 2
+   }
+   
+   sigma <- noise_sd_base * pcr_duplicates / sqrt(coverage / 10)
+   observed_cr <- true_cr + rnorm(1, 0, sigma)
+   log2_ratio <- log2(observed_cr)
+   return(log2_ratio)
+ }
> # ============================================
> # SCENARIO B SIMULATION (custom noise)
> # ============================================
> 
> simulate_scenario_B_custom <- function(true_mosaicism, n_chr = 24) {
+   affected_idx <- sample(1:n_chr, 1)
+   true_log2 <- rep(0, n_chr)
+   true_log2[affected_idx] <- log2(1 + true_mosaicism / 2)
+   
+   sigma <- noise_sd_base * pcr_duplicates / sqrt(coverage / 10)
+   chr_noise <- rnorm(n_chr, 0, sigma)
+   batch_effect <- if (correlated_noise) rnorm(1, 0, batch_sd) else 0
+   
+   observed <- true_log2 + chr_noise + batch_effect
+   return(observed)
+ }
> # ============================================
> # RUN SINGLE CHROMOSOME (if selected)
> # ============================================
> 
> if (run_single_chr) {
+   cat("\n===== Single chromosome =====\n")
+   
+   # M = 0% (false positive rate)
+   vals0 <- replicate(n_replicates, simulate_cnv_custom(0))
+   fp_rate <- mean(abs(vals0) >= low_thr)
+   
+   # M = 15% (false negative rate = 1 - sensitivity)
+   vals15 <- replicate(n_replicates, simulate_cnv_custom(0.15))
+   fn_rate <- mean(abs(vals15) < low_thr)
+   
+   cat("False positive rate (M=0%):", round(fp_rate, 4), "\n")
+   cat("False negative rate (M=15%):", round(fn_rate, 4), "\n")
+ }

===== Single chromosome =====
False positive rate (M=0%): 0.086 
False negative rate (M=15%): 0.724 
> # ============================================
> # RUN SCENARIO B (if selected)
> # ============================================
> 
> if (run_scenario_B) {
+   cat("\n===== Scenario B (single affected chromosome, 24 total) =====\n")
+   
+   # M = 0% (false positive rate)
+   fp_count <- 0
+   for (i in 1:n_replicates) {
+     obs <- simulate_scenario_B_custom(0)
+     classes <- classify_embryo(obs, low_thr, high_thr)
+     if (any(classes != "normal")) fp_count <- fp_count + 1
+   }
+   fp_rate_B <- fp_count / n_replicates
+   
+   # M = 15% (false negative rate)
+   fn_count <- 0
+   for (i in 1:n_replicates) {
+     obs <- simulate_scenario_B_custom(0.15)
+     classes <- classify_embryo(obs, low_thr, high_thr)
+     if (all(classes == "normal")) fn_count <- fn_count + 1
+   }
+   fn_rate_B <- fn_count / n_replicates
+   
+   cat("False positive rate (M=0%):", round(fp_rate_B, 4), "\n")
+   cat("False negative rate (M=15%):", round(fn_rate_B, 4), "\n")
+ }

===== Scenario B (single affected chromosome, 24 total) =====
False positive rate (M=0%): 0.539 
False negative rate (M=15%): 0.348 
> 
