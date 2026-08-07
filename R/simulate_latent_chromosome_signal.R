# Simulate latent chromosome-level signal for the realistic aggregated-signal
# model.

simulate_latent_chromosome_signal <- function(
    embryo_state,
    chromosomes = paste0("chr", 1:22),
    normal_spurious_abnormal_prob = 0.10,
    ambiguous_n_abnormal = c(1L, 2L),
    abnormal_n_abnormal = 2:6,
    ambiguous_effect_mean = 0.22,
    ambiguous_effect_sd = 0.05,
    abnormal_effect_mean = 0.48,
    abnormal_effect_sd = 0.10) {
  out <- vector("list", nrow(embryo_state))

  for (i in seq_len(nrow(embryo_state))) {
    status_i <- embryo_state$latent_status[i]
    n_chr <- length(chromosomes)

    latent_signal <- rep(0, n_chr)
    chromosome_state <- rep("baseline", n_chr)

    if (status_i == "normal") {
      n_affected <- rbinom(1, size = 1, prob = normal_spurious_abnormal_prob)
      effect_mean <- ambiguous_effect_mean * 0.75
      effect_sd <- ambiguous_effect_sd
    } else if (status_i == "ambiguous") {
      n_affected <- sample(ambiguous_n_abnormal, size = 1)
      effect_mean <- ambiguous_effect_mean
      effect_sd <- ambiguous_effect_sd
    } else if (status_i == "abnormal") {
      n_affected <- sample(abnormal_n_abnormal, size = 1)
      effect_mean <- abnormal_effect_mean
      effect_sd <- abnormal_effect_sd
    } else {
      stop("Unknown latent_status in embryo_state.", call. = FALSE)
    }

    if (n_affected > 0) {
      idx <- sample(seq_len(n_chr), size = n_affected, replace = FALSE)
      effects <- rnorm(n_affected, mean = effect_mean, sd = effect_sd)
      effects <- pmax(effects, 0.05)
      signs <- sample(c(-1, 1), size = n_affected, replace = TRUE)

      latent_signal[idx] <- effects * signs
      chromosome_state[idx] <- if (status_i == "ambiguous") "ambiguous_event" else "abnormal_event"
    }

    out[[i]] <- data.frame(
      embryo_id = embryo_state$embryo_id[i],
      latent_status = status_i,
      chromosome = chromosomes,
      chromosome_state = chromosome_state,
      latent_signal = latent_signal,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}
