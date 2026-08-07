
run_scenario_A <- function(m_levels = c(0, 0.15, 0.20),
                           n_chr = 24, coverage = 10,
                           low_thr = 0.15, high_thr = 0.30,
                           n_replicates = 5000) {
  
  source("R/simulate_embryo_uniform.R")
  source("R/classify_embryo.R")
  
  results <- data.frame()
  
  for (m in m_levels) {
    for (rep in 1:n_replicates) {
      obs <- simulate_embryo_uniform(m, n_chr = n_chr, coverage = coverage)
      classes <- classify_embryo(obs, low_thr, high_thr)
      embryo_abnormal <- any(classes != "normal")
      results <- rbind(results, data.frame(m = m, abnormal = embryo_abnormal))
    }
  }
  
  agg <- aggregate(abnormal ~ m, data = results, FUN = mean)
  colnames(agg)[2] <- "prop_abnormal"
  agg$prop_normal <- 1 - agg$prop_abnormal
  
  return(agg)
}

