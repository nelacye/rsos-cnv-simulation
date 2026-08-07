
run_scenario_single <- function(m_levels = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30),
                                coverages = c(5, 10, 20, 50),
                                low_thr = 0.15, high_thr = 0.30,
                                n_replicates = 5000) {
  
  source("R/simulate_cnv.R")
  source("R/classify_embryo.R")
  
  results <- expand.grid(m = m_levels, cov = coverages, rep = 1:n_replicates)
  results$log2r <- NA
  
  for (i in 1:nrow(results)) {
    results$log2r[i] <- simulate_cnv(results$m[i], coverage = results$cov[i])
  }
  
  results$class <- classify_embryo(results$log2r, low_thr, high_thr)
  
  agg <- aggregate(class ~ m + cov, data = results,
                   FUN = function(x) mean(x == "normal"))
  colnames(agg)[3] <- "prop_normal"
  
  return(agg)
}

