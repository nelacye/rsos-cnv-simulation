# Comparable idealized simulation that preserves the original idealized signal
# model while allowing the latent chromosome burden to match the realistic
# branch.

simulate_embryo_comparable_idealized <- function(
    true_mosaicism,
    n_affected = 1L,
    n_chr = 24,
    coverage = 10,
    noise_sd_base = 0.05,
    pcr_duplicates = 1.2,
    batch_sd = 0.05,
    correlated_noise = TRUE) {
  stopifnot(n_affected >= 0, n_affected <= n_chr)

  true_log2 <- rep(0, n_chr)

  if (n_affected > 0) {
    affected_idx <- sample(seq_len(n_chr), size = n_affected, replace = FALSE)
    true_log2[affected_idx] <- log2(1 + true_mosaicism / 2)
  }

  chr_noise <- rnorm(n_chr, 0, noise_sd_base * pcr_duplicates / sqrt(coverage / 10))
  batch_effect <- if (correlated_noise) rnorm(1, 0, batch_sd) else 0

  true_log2 + chr_noise + batch_effect
}
