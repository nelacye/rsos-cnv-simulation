# Classification helpers for the realistic aggregated-signal model.

estimate_empirical_signal_cutoff <- function(observed_data, quantile_cutoff = 0.95) {
  quantile(abs(observed_data$observed_signal), probs = quantile_cutoff, na.rm = TRUE)
}

posterior_abnormal_signal <- function(
    observed_signal,
    coverage,
    prior_normal = 0.85,
    prior_ambiguous = 0.10,
    prior_abnormal = 0.05,
    reference_coverage = 10,
    base_sd = 0.28,
    ambiguous_effect = 0.18,
    abnormal_effect = 0.38,
    component_sd = 0.15) {
  effective_sd <- base_sd * sqrt(reference_coverage / coverage)
  ambiguous_sd <- sqrt(effective_sd^2 + component_sd^2)
  abnormal_sd <- sqrt(effective_sd^2 + component_sd^2)

  likelihood_normal <- dnorm(observed_signal, mean = 0, sd = effective_sd)
  likelihood_ambiguous <- 0.5 * dnorm(observed_signal, mean = ambiguous_effect, sd = ambiguous_sd) +
    0.5 * dnorm(observed_signal, mean = -ambiguous_effect, sd = ambiguous_sd)
  likelihood_abnormal <- 0.5 * dnorm(observed_signal, mean = abnormal_effect, sd = abnormal_sd) +
    0.5 * dnorm(observed_signal, mean = -abnormal_effect, sd = abnormal_sd)

  denominator <- prior_normal * likelihood_normal +
    prior_ambiguous * likelihood_ambiguous +
    prior_abnormal * likelihood_abnormal

  numerator <- prior_ambiguous * likelihood_ambiguous +
    prior_abnormal * likelihood_abnormal

  pmin(pmax(numerator / denominator, 0), 1)
}

compute_chromosome_calls_realistic <- function(
    observed_data,
    chromosome_cutoff,
    posterior_cutoff = 0.50) {
  out <- observed_data

  out$threshold_chromosome_call <- ifelse(
    abs(out$observed_signal) > chromosome_cutoff,
    "abnormal",
    "normal"
  )

  out$posterior_abnormal <- posterior_abnormal_signal(
    observed_signal = out$observed_signal,
    coverage = out$coverage,
    reference_coverage = unique(out$reference_coverage)[1]
  )

  out$bayes_chromosome_call <- ifelse(
    out$posterior_abnormal > posterior_cutoff,
    "abnormal",
    "normal"
  )

  out
}

summarize_embryo_realistic <- function(chromosome_calls) {
  embryo_ids <- unique(chromosome_calls$embryo_id)
  out <- vector("list", length(embryo_ids))

  for (i in seq_along(embryo_ids)) {
    emb_id <- embryo_ids[i]
    emb <- chromosome_calls[chromosome_calls$embryo_id == emb_id, ]

    out[[i]] <- data.frame(
      embryo_id = emb_id,
      latent_status = emb$latent_status[1],
      coverage = emb$coverage[1],
      n_chromosomes = nrow(emb),
      threshold_n_abnormal = sum(emb$threshold_chromosome_call == "abnormal"),
      bayes_n_abnormal = sum(emb$bayes_chromosome_call == "abnormal"),
      mean_abs_signal = mean(abs(emb$observed_signal)),
      sum_abs_signal = sum(abs(emb$observed_signal)),
      mean_posterior = mean(emb$posterior_abnormal),
      sum_posterior = sum(emb$posterior_abnormal),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

estimate_aggregate_cutoffs_realistic <- function(
    calibration_summary,
    normal_status = "normal",
    cutoff_quantile = 0.95,
    evaluation_label = "out-of-sample evaluation") {
  normal_summary <- calibration_summary[calibration_summary$latent_status == normal_status, ]

  if (nrow(normal_summary) == 0) {
    stop("Calibration summary must contain at least one normal embryo.", call. = FALSE)
  }

  data.frame(
    evaluation_label = evaluation_label,
    calibration_n_embryos = nrow(calibration_summary),
    calibration_n_normal = nrow(normal_summary),
    threshold_mean_cutoff = quantile(
      normal_summary$mean_abs_signal,
      probs = cutoff_quantile,
      na.rm = TRUE
    ),
    threshold_sum_cutoff = quantile(
      normal_summary$sum_abs_signal,
      probs = cutoff_quantile,
      na.rm = TRUE
    ),
    bayes_mean_cutoff = quantile(
      normal_summary$mean_posterior,
      probs = cutoff_quantile,
      na.rm = TRUE
    ),
    bayes_sum_cutoff = quantile(
      normal_summary$sum_posterior,
      probs = cutoff_quantile,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

apply_embryo_decision_rules_realistic <- function(
    embryo_summary,
    aggregate_cutoffs = NULL,
    evaluation_label = "out-of-sample evaluation") {
  rules <- rbind(
    data.frame(
      classifier = "threshold",
      rule_name = c("ge_1_abnormal_chromosome", "ge_2_abnormal_chromosomes", "ge_3_abnormal_chromosomes"),
      embryo_call = I(list(
        ifelse(embryo_summary$threshold_n_abnormal >= 1, "abnormal", "normal"),
        ifelse(embryo_summary$threshold_n_abnormal >= 2, "abnormal", "normal"),
        ifelse(embryo_summary$threshold_n_abnormal >= 3, "abnormal", "normal")
      )),
      stringsAsFactors = FALSE
    ),
    data.frame(
      classifier = "bayesian",
      rule_name = c("ge_1_abnormal_chromosome", "ge_2_abnormal_chromosomes", "ge_3_abnormal_chromosomes"),
      embryo_call = I(list(
        ifelse(embryo_summary$bayes_n_abnormal >= 1, "abnormal", "normal"),
        ifelse(embryo_summary$bayes_n_abnormal >= 2, "abnormal", "normal"),
        ifelse(embryo_summary$bayes_n_abnormal >= 3, "abnormal", "normal")
      )),
      stringsAsFactors = FALSE
    )
  )

  if (is.null(aggregate_cutoffs)) {
    stop("aggregate_cutoffs must be supplied for realistic embryo aggregation rules.", call. = FALSE)
  }

  threshold_mean_cutoff <- aggregate_cutoffs$threshold_mean_cutoff[1]
  threshold_sum_cutoff <- aggregate_cutoffs$threshold_sum_cutoff[1]
  bayes_mean_cutoff <- aggregate_cutoffs$bayes_mean_cutoff[1]
  bayes_sum_cutoff <- aggregate_cutoffs$bayes_sum_cutoff[1]

  rules <- rbind(
    rules,
    data.frame(
      classifier = "threshold",
      rule_name = c("mean_absolute_signal", "sum_absolute_signal"),
      embryo_call = I(list(
        ifelse(embryo_summary$mean_abs_signal > threshold_mean_cutoff, "abnormal", "normal"),
        ifelse(embryo_summary$sum_abs_signal > threshold_sum_cutoff, "abnormal", "normal")
      )),
      stringsAsFactors = FALSE
    ),
    data.frame(
      classifier = "bayesian",
      rule_name = c("mean_posterior_signal", "sum_posterior_signal"),
      embryo_call = I(list(
        ifelse(embryo_summary$mean_posterior > bayes_mean_cutoff, "abnormal", "normal"),
        ifelse(embryo_summary$sum_posterior > bayes_sum_cutoff, "abnormal", "normal")
      )),
      stringsAsFactors = FALSE
    )
  )

  out <- vector("list", nrow(rules))

  for (i in seq_len(nrow(rules))) {
    tmp <- embryo_summary
    tmp$evaluation_label <- evaluation_label
    tmp$classifier <- rules$classifier[i]
    tmp$rule_name <- rules$rule_name[i]
    tmp$embryo_call <- rules$embryo_call[[i]]
    out[[i]] <- tmp
  }

  list(
    embryo_calls = do.call(rbind, out),
    aggregate_thresholds = aggregate_cutoffs
  )
}

score_realistic_decisions <- function(embryo_calls) {
  scored <- embryo_calls
  scored$truth_binary <- ifelse(scored$latent_status == "normal", "normal", "abnormal")

  split_keys <- unique(scored[, c("classifier", "rule_name")])
  metrics_out <- vector("list", nrow(split_keys))
  confusion_out <- vector("list", nrow(split_keys))

  for (i in seq_len(nrow(split_keys))) {
    classifier_i <- split_keys$classifier[i]
    rule_i <- split_keys$rule_name[i]
    tmp <- scored[scored$classifier == classifier_i & scored$rule_name == rule_i, ]

    tp <- sum(tmp$truth_binary == "abnormal" & tmp$embryo_call == "abnormal")
    tn <- sum(tmp$truth_binary == "normal" & tmp$embryo_call == "normal")
    fp <- sum(tmp$truth_binary == "normal" & tmp$embryo_call == "abnormal")
    fn <- sum(tmp$truth_binary == "abnormal" & tmp$embryo_call == "normal")

    metrics_out[[i]] <- data.frame(
      evaluation_label = tmp$evaluation_label[1],
      classifier = classifier_i,
      rule_name = rule_i,
      sensitivity = tp / (tp + fn),
      specificity = tn / (tn + fp),
      accuracy = (tp + tn) / nrow(tmp),
      true_positive = tp,
      true_negative = tn,
      false_positive = fp,
      false_negative = fn,
      stringsAsFactors = FALSE
    )

    confusion_out[[i]] <- aggregate(
      embryo_id ~ truth_binary + embryo_call,
      data = tmp,
      FUN = length
    )
    names(confusion_out[[i]])[3] <- "n"
    confusion_out[[i]]$evaluation_label <- tmp$evaluation_label[1]
    confusion_out[[i]]$classifier <- classifier_i
    confusion_out[[i]]$rule_name <- rule_i
  }

  list(
    metrics = do.call(rbind, metrics_out),
    confusion = do.call(rbind, confusion_out)
  )
}
