# Full numerical reconstruction. Allow time for all analysis commands.
# Run only when ready to perform the clean-session validation described in README.
# Existing project data, objects and outputs are preserved.
local({
  recovery_dir <- "recovered_analysis"
  run_mode <- "Reconstructed numeric pipeline plus historical sensitivity checks when available"
  source(file.path(recovery_dir, "R/00_setup.R"), local = TRUE)
  source(file.path(recovery_dir, "R/01_reconstructed_baseline.R"), local = TRUE)
  if (file.exists(checkpoint) || file.exists("reference/historical_cell_labels.csv")) {
    pbmc_annotated <- load_historical_baseline()
    source(file.path(recovery_dir, "R/02_marker_review.R"), local = TRUE)
    source(file.path(recovery_dir, "R/03_sensitivity_checks.R"), local = TRUE)
    source(file.path(recovery_dir, "R/04_saved_baseline_figures.R"), local = TRUE)
  } else {
    message("Numeric reconstruction finished; historical checkpoint absent. ",
            "Historical annotation/sensitivity/figure restoration was skipped. ",
            "Review markers before adding candidate labels to rebuilt clusters.")
  }
  source(file.path(recovery_dir, "R/05_environment_record.R"), local = TRUE)
})
