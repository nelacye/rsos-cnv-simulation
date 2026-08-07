
# classify_embryo.R
# Classify one chromosome based on log2 ratio and thresholds

classify_embryo <- function(log2r, low_thr, high_thr) {
  ifelse(abs(log2r) < low_thr, "normal",
         ifelse(abs(log2r) < high_thr, "mosaic", "aneuploid"))
}

