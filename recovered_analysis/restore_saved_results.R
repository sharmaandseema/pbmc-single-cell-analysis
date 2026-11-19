# Lightweight restoration: uses your existing saved annotated object.
# No normalization, PCA, clustering or UMAP is recomputed.
# Run from PBMC_project with:
# source("recovered_analysis/restore_saved_results.R")
local({
  recovery_dir <- "recovered_analysis"
  run_mode <- "Restored figures/tables from historical annotated checkpoint"
  source(file.path(recovery_dir, "R/00_setup.R"), local = TRUE)
  pbmc_annotated <- load_historical_baseline()
  source(file.path(recovery_dir, "R/02_marker_review.R"), local = TRUE)
  source(file.path(recovery_dir, "R/04_saved_baseline_figures.R"), local = TRUE)
  source(file.path(recovery_dir, "R/05_environment_record.R"), local = TRUE)
})
