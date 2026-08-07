# Convert latent chromosome signal into an observed pseudo log2-like signal.
#
# The goal is not to reconstruct QDNAseq exactly. The goal is to introduce an
# explicit distortion / aggregation layer that can reproduce the empirical
# pattern we saw in workbook-derived signal:
#   - autosomes roughly centered after within-embryo normalization
#   - moderate sensitivity
#   - weak specificity under simple embryo-level rules

simulate_measurement_layer <- function(
    latent_chromosome_signal,
    embryo_state,
    raw_baseline = 1.70,
    base_noise_sd = 0.16,
    chromosome_scale_sdlog = 0.18,
    status_noise_multiplier = c(normal = 1.15, ambiguous = 1.00, abnormal = 0.95),
    min_raw_signal = 1e-4) {
  dat <- merge(
    latent_chromosome_signal,
    embryo_state,
    by = c("embryo_id", "latent_status"),
    all.x = TRUE,
    sort = FALSE
  )

  dat$chromosome_scale <- rlnorm(nrow(dat), meanlog = 0, sdlog = chromosome_scale_sdlog)
  dat$coverage_noise_sd <- base_noise_sd * sqrt(dat$reference_coverage / dat$coverage)

  dat$status_noise_multiplier <- status_noise_multiplier[dat$latent_status]
  dat$raw_signal_mean <- raw_baseline * dat$embryo_scale * dat$chromosome_scale * (2 ^ dat$latent_signal)
  dat$raw_signal <- dat$raw_signal_mean + dat$baseline_shift +
    rnorm(nrow(dat), mean = 0, sd = dat$coverage_noise_sd * dat$status_noise_multiplier)
  dat$raw_signal <- pmax(dat$raw_signal, min_raw_signal)

  autosome_median <- aggregate(
    raw_signal ~ embryo_id,
    data = dat,
    FUN = median
  )
  names(autosome_median)[2] <- "autosome_median_raw"

  dat <- merge(dat, autosome_median, by = "embryo_id", all.x = TRUE, sort = FALSE)
  dat$normalized_signal <- dat$raw_signal / dat$autosome_median_raw
  dat$observed_signal <- log2(dat$normalized_signal)
  dat$abs_observed_signal <- abs(dat$observed_signal)

  dat[order(dat$embryo_id, dat$chromosome), ]
}
