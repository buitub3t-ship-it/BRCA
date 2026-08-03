source("config.R")
source(file.path("R", "author_plots.R"))

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required. Run 00_install_packages.R first.")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required. Run 00_install_packages.R first.")
}

required_files <- c(
  internal_metrics = file.path(RESULTS_DIR, "03_internal_full_monte_carlo_metrics.csv"),
  external_metrics = file.path(RESULTS_DIR, "04_expression_external_monte_carlo_metrics.csv"),
  final_external = file.path(RESULTS_DIR, "04_final_expression_external_model.rds")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Required training results are missing:\n",
    paste("-", missing_files, collapse = "\n"),
    "\nRun scripts 03 and 04 first."
  )
}

internal_metrics <- data.table::fread(required_files[["internal_metrics"]], data.table = FALSE)
external_metrics <- data.table::fread(required_files[["external_metrics"]], data.table = FALSE)
final_external_object <- readRDS(required_files[["final_external"]])
final_external <- final_external_object$result

# Figure 2A-like: full TCGA Meth+Exp internal validation.
for (ext in c("pdf", "png")) {
  save_author_figure2a(
    internal_metrics,
    file.path(RESULTS_DIR, paste0("06_Figure2A_internal_Cindex.", ext))
  )
}

# Figure 2B-like: expression-only internal/external validation.
for (ext in c("pdf", "png")) {
  save_author_figure2b(
    external_metrics,
    file.path(RESULTS_DIR, paste0("06_Figure2B_external_Cindex.", ext))
  )
}

# Figure 3-like: 10-year Kaplan-Meier curves for quartile-defined risk groups.
for (ext in c("pdf", "png")) {
  grouped_predictions <- save_author_figure3(
    final_external$predictions,
    file.path(RESULTS_DIR, paste0("06_Figure3_KaplanMeier_risk_groups.", ext))
  )
}

risk_assignment <- do.call(
  rbind,
  lapply(names(grouped_predictions), function(cohort) {
    d <- grouped_predictions[[cohort]]
    d$cohort <- sub("_train$", "", cohort)
    d
  })
)
rownames(risk_assignment) <- NULL
data.table::fwrite(
  risk_assignment,
  file.path(RESULTS_DIR, "06_risk_group_assignments.csv")
)

hr_table <- do.call(
  rbind,
  lapply(names(grouped_predictions), function(cohort) {
    risk_group_hr_table(
      grouped_predictions[[cohort]],
      sub("_train$", "", cohort)
    )
  })
)
data.table::fwrite(
  hr_table,
  file.path(RESULTS_DIR, "06_risk_group_hazard_ratios.csv")
)

# Figure S1A-like: grouped 5-year calibration of the final M2EFM Exp+Clin model.
calibration_table <- build_calibration_table(
  final_model = final_external$fit$final_model,
  predictions = final_external$predictions,
  horizon_days = 5 * 365.25,
  n_bins = 5L
)
data.table::fwrite(
  calibration_table,
  file.path(RESULTS_DIR, "06_calibration_5year.csv")
)
for (ext in c("pdf", "png")) {
  save_author_figureS1(
    calibration_table,
    file.path(RESULTS_DIR, paste0("06_FigureS1_calibration_5year.", ext))
  )
}

message("Author-style plots completed. See ./results files beginning with 06_.")
