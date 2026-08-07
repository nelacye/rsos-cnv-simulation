
# simulate_cnv.R
# Simulate observed log2 ratio for one CNV on one chromosome

simulate_cnv <- function(true_mosaicism, cnv_type = "dup",
                         coverage = 10, noise_sd_base = 0.05,
                         pcr_duplicates = 1.2) {
  
  if (cnv_type == "dup") {
    true_cr <- 1 + true_mosaicism / 2
  } else if (cnv_type == "del") {
    true_cr <- 1 - true_mosaicism / 2
  } else {
    stop("cnv_type must be 'dup' or 'del'")
  }
  
  effective_noise_sd <- noise_sd_base * pcr_duplicates / sqrt(coverage / 10)
  observed_cr <- true_cr + rnorm(1, 0, effective_noise_sd)
  log2_ratio <- log2(observed_cr)
  
  return(log2_ratio)
}

