
# simulate_embryo_single_affected.R
# Scenario B: only one random chromosome carries the mosaic CNV

simulate_embryo_single_affected <- function(true_mosaicism, n_chr = 24, coverage = 10,
                                            noise_sd_base = 0.05, pcr_duplicates = 1.2,
                                            batch_sd = 0.05, correlated_noise = TRUE) {
  
  affected_idx <- sample(1:n_chr, 1)
  true_log2 <- rep(0, n_chr)
  true_log2[affected_idx] <- log2(1 + true_mosaicism / 2)
  chr_noise <- rnorm(n_chr, 0, noise_sd_base * pcr_duplicates / sqrt(coverage / 10))
  batch_effect <- if (correlated_noise) rnorm(1, 0, batch_sd) else 0
  observed <- true_log2 + chr_noise + batch_effect
  
  return(observed)
}

