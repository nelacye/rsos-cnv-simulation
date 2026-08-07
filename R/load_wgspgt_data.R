# Helpers for first-pass WGS-PGT validation data preparation.
#
# This loader intentionally reads only the two workbook sheets needed for the
# current descriptive validation workflow:
#   - Fig1f_WGS for chromosome-level measurements
#   - Fig1d_Fig1e for embryo-level metadata

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
})

coerce_decimal_numeric <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  x_chr <- as.character(x)
  x_chr <- stringr::str_trim(x_chr)
  x_chr[x_chr %in% c("", "NA", "N/A", "na", "null", "NULL")] <- NA_character_
  x_chr <- stringr::str_replace_all(x_chr, "\\s+", "")

  both_marks <- stringr::str_detect(x_chr, ",") & stringr::str_detect(x_chr, "\\.")
  comma_decimal <- both_marks & stringr::str_detect(x_chr, ",\\d+$")
  dot_decimal <- both_marks & stringr::str_detect(x_chr, "\\.\\d+$")

  x_chr[comma_decimal] <- stringr::str_replace_all(x_chr[comma_decimal], "\\.", "")
  x_chr[comma_decimal] <- stringr::str_replace_all(x_chr[comma_decimal], ",", ".")
  x_chr[dot_decimal] <- stringr::str_replace_all(x_chr[dot_decimal], ",", "")

  comma_only <- !both_marks & stringr::str_detect(x_chr, ",")
  x_chr[comma_only] <- stringr::str_replace_all(x_chr[comma_only], ",", ".")

  suppressWarnings(as.numeric(x_chr))
}

load_wgspgt_sheet_data <- function(
    workbook_path = file.path("data", "wgspgt", "Fig.1d_Fig.1e_Fig.1f_Fig.1g_Fig.2b_Fig5.xlsx"),
    wgs_sheet = "Fig1f_WGS",
    metadata_sheet = "Fig1d_Fig1e") {
  if (!file.exists(workbook_path)) {
    stop(sprintf("Workbook does not exist: %s", workbook_path), call. = FALSE)
  }

  wgs_raw <- readxl::read_excel(workbook_path, sheet = wgs_sheet, .name_repair = "unique_quiet")
  metadata_raw <- readxl::read_excel(workbook_path, sheet = metadata_sheet, .name_repair = "unique_quiet")

  wgs_tidy <- wgs_raw %>%
    mutate(EmbryoNumber = as.character(EmbryoNumber)) %>%
    pivot_longer(
      cols = matches("^chr"),
      names_to = "chromosome",
      values_to = "log2r"
    ) %>%
    transmute(
      embryo_id = EmbryoNumber,
      chromosome = as.character(chromosome),
      log2r = coerce_decimal_numeric(log2r),
      abs_log2r = abs(log2r)
    )

  metadata_tidy <- metadata_raw %>%
    filter(!is.na(EmbryoNumber), EmbryoNumber != "") %>%
    transmute(
      embryo_id = as.character(EmbryoNumber),
      meanDepth = coerce_decimal_numeric(meanDepth),
      meanBreadth = coerce_decimal_numeric(meanBreadth),
      EmbryoStatus = as.character(EmbryoStatus),
      Method = as.character(Method)
    ) %>%
    distinct(embryo_id, .keep_all = TRUE)

  wgs_tidy %>%
    left_join(metadata_tidy, by = "embryo_id") %>%
    arrange(embryo_id, chromosome)
}

save_wgspgt_cleaned_data <- function(
    workbook_path = file.path("data", "wgspgt", "Fig.1d_Fig.1e_Fig.1f_Fig.1g_Fig.2b_Fig5.xlsx"),
    output_path = file.path("results", "wgspgt_validation", "wgspgt_cleaned.csv")) {
  cleaned <- load_wgspgt_sheet_data(workbook_path = workbook_path)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(cleaned, output_path)

  cleaned
}

read_wgspgt_cleaned_data <- function(
    cleaned_path = file.path("results", "wgspgt_validation", "wgspgt_cleaned.csv")) {
  if (!file.exists(cleaned_path)) {
    stop(sprintf("Cleaned validation file does not exist: %s", cleaned_path), call. = FALSE)
  }

  readr::read_csv(cleaned_path, show_col_types = FALSE) %>%
    mutate(
      embryo_id = as.character(embryo_id),
      chromosome = as.character(chromosome),
      EmbryoStatus = as.character(EmbryoStatus),
      Method = as.character(Method)
    )
}
