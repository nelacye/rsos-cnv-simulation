
# simulate_embryo_uniform.R
# Scenario A: all chromosomes have the same true mosaicism level

simulate_embryo_uniform <- function(true_mosaicism, n_chr = 24, coverage = 10,
                                    noise_sd_base = 0.05, pcr_duplicates = 1.2,
                                    batch_sd = 0.05, correlated_noise = TRUE) {
  
  true_log2 <- log2(1 + true_mosaicism / 2)
  chr_noise <- rnorm(n_chr, 0, noise_sd_base * pcr_duplicates / sqrt(coverage / 10))
  batch_effect <- if (correlated_noise) rnorm(1, 0, batch_sd) else 0
  observed <- rep(true_log2, n_chr) + chr_noise + batch_effect
  
  return(observed)
}

