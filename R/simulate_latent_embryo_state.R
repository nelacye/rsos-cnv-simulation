# Simulate embryo-level latent states for the realistic aggregated-signal model.

simulate_latent_embryo_state <- function(
    n_embryos,
    status_probs = c(normal = 0.35, ambiguous = 0.30, abnormal = 0.35),
    coverage_meanlog = log(10),
    coverage_sdlog = 0.35,
    baseline_shift_sd = 0.10,
    embryo_scale_sdlog = 0.18,
    reference_coverage = 10) {
  stopifnot(abs(sum(status_probs) - 1) < 1e-8)

  latent_status <- sample(
    x = names(status_probs),
    size = n_embryos,
    replace = TRUE,
    prob = status_probs
  )

  coverage <- rlnorm(n_embryos, meanlog = coverage_meanlog, sdlog = coverage_sdlog)
  coverage <- pmax(3, pmin(40, coverage))

  data.frame(
    embryo_id = sprintf("sim_%04d", seq_len(n_embryos)),
    latent_status = latent_status,
    coverage = coverage,
    reference_coverage = reference_coverage,
    baseline_shift = rnorm(n_embryos, mean = 0, sd = baseline_shift_sd),
    embryo_scale = rlnorm(n_embryos, meanlog = 0, sdlog = embryo_scale_sdlog),
    stringsAsFactors = FALSE
  )
}
