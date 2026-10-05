# Shared data prep for the REPO1 compound figures. Paths are relative to a notebook in notebooks/.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
})
source(file.path("..", "plotting_helpers", "repo1_themes.R"))

plate_map_file <- file.path("..", "REPO1_April2025.xlsx")

# One row per unique compound: wells without a compound are dropped and compounds plated in two wells are counted once.
# Phase is grouped into the levels in phase_levels (the clinical trial phases are combined; missing becomes "Not annotated").
read_repo1_compounds <- function(path = plate_map_file) {
  read_excel(path, sheet = "REPO2025_PlateMap") %>%
    filter(!is.na(BROAD_CPD_ID)) %>%
    distinct(BROAD_CPD_ID, .keep_all = TRUE) %>%
    mutate(phase = factor(case_when(
      is.na(Phase) ~ "Not annotated",
      str_detect(Phase, "^Phase") ~ "Phase 1-3",
      TRUE ~ Phase
    ), levels = phase_levels))
}

# Split a multi-valued annotation column (comma separated) into one row per compound and value.
# Values are trimmed and lower-cased so spelling variants such as "Oncology" and "oncology" are counted together.
explode_annotation <- function(compounds, column) {
  compounds %>%
    filter(!is.na(.data[[column]])) %>%
    mutate(value = str_split(.data[[column]], ",")) %>%
    unnest(value) %>%
    mutate(value = str_to_lower(str_squish(value))) %>%
    filter(value != "") %>%
    distinct(BROAD_CPD_ID, value, .keep_all = TRUE)
}
