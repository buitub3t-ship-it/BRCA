source("config.R")
source(file.path("R", "io_and_preprocessing.R"))
source(file.path("R", "modeling.R"))

needed <- c(
  "tcga_exp_common.rds", "kao_exp_common.rds", "teru_exp_common.rds",
  "tcga_clin_clean.rds", "kao_clin_clean.rds", "teru_clin_clean.rds",
  "m2eqtl_discovery.rds"
)
if (!all(file.exists(file.path(CACHE_DIR, needed)))) {
  stop("Run scripts 01 and 02 first.")
}

tcga_exp <- readRDS(file.path(CACHE_DIR, "tcga_exp_common.rds"))
kao_exp <- readRDS(file.path(CACHE_DIR, "kao_exp_common.rds"))
teru_exp <- readRDS(file.path(CACHE_DIR, "teru_exp_common.rds"))
tcga_clin <- readRDS(file.path(CACHE_DIR, "tcga_clin_clean.rds"))
kao_clin <- readRDS(file.path(CACHE_DIR, "kao_clin_clean.rds"))
teru_clin <- readRDS(file.path(CACHE_DIR, "teru_clin_clean.rds"))
m2e <- readRDS(file.path(CACHE_DIR, "m2eqtl_discovery.rds"))

signature <- m2e$expression_signature
missing_by_cohort <- list(
  TCGA = setdiff(signature, rownames(tcga_exp)),
  Kao = setdiff(signature, rownames(kao_exp)),
  Terunuma = setdiff(signature, rownames(teru_exp))
)
if (any(lengths(missing_by_cohort) > 0L)) {
  stop("Expression signature is not available in every cohort.")
}

tcga_sig <- tcga_exp[signature, , drop = FALSE]
kao_sig <- kao_exp[signature, , drop = FALSE]
teru_sig <- teru_exp[signature, , drop = FALSE]

# Match the author's evaluate(): learn mean/range from TCGA and apply to each
# external validation cohort.
tcga_scale <- scale_features_reference(tcga_sig)
tcga_scaled <- tcga_scale$values
kao_scaled <- apply_feature_reference(
  kao_sig,
  tcga_scale$means,
  tcga_scale$ranges
)
teru_scaled <- apply_feature_reference(
  teru_sig,
  tcga_scale$means,
  tcga_scale$ranges
)

validation_sets <- list(
  Kao = list(features = kao_scaled, clinical = kao_clin),
  Terunuma = list(features = teru_scaled, clinical = teru_clin)
)

external <- monte_carlo_m2efm(
  training_feature_by_sample = tcga_scaled,
  training_clinical = tcga_clin,
  validation_sets = validation_sets,
  n_splits = N_MONTE_CARLO_SPLITS
)

saveRDS(external, file.path(CACHE_DIR, "external_expression_monte_carlo.rds"))
data.table::fwrite(
  external$metrics,
  file.path(RESULTS_DIR, "04_expression_external_monte_carlo_metrics.csv")
)

summary_metrics <- summarize_cindex_metrics(external$metrics)
data.table::fwrite(
  summary_metrics,
  file.path(RESULTS_DIR, "04_expression_external_summary.csv")
)

save_metric_boxplot(
  external$metrics,
  file.path(RESULTS_DIR, "04_expression_external_boxplot.pdf"),
  "M2EFM Exp+Clin: TCGA internal and external validation"
)

# Final expression-only model trained on all eligible TCGA samples.
final_external <- fit_final_and_validate(
  training_feature_by_sample = tcga_scaled,
  training_clinical = tcga_clin,
  validation_sets = validation_sets,
  seed = 20001L
)

saveRDS(
  list(
    result = final_external,
    scaling = tcga_scale[c("means", "ranges")],
    expression_signature = signature
  ),
  file.path(RESULTS_DIR, "04_final_expression_external_model.rds")
)
data.table::fwrite(
  final_external$metrics,
  file.path(RESULTS_DIR, "04_final_expression_external_metrics.csv")
)

for (cohort in names(final_external$predictions)) {
  data.table::fwrite(
    final_external$predictions[[cohort]],
    file.path(RESULTS_DIR, paste0("04_predictions_", cohort, ".csv")),
    row.names = FALSE
  )
}

message("Expression-only external-validation training completed.")
print(summary_metrics)
print(final_external$metrics)
