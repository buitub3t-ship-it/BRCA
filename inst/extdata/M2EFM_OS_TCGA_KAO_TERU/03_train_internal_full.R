source("config.R")
source(file.path("R", "io_and_preprocessing.R"))
source(file.path("R", "m2eqtl.R"))
source(file.path("R", "modeling.R"))

needed <- c(
  "tcga_exp_common.rds", "tcga_meth_beta.rds",
  "tcga_clin_clean.rds", "m2eqtl_discovery.rds"
)
if (!all(file.exists(file.path(CACHE_DIR, needed)))) {
  stop("Run scripts 01 and 02 first.")
}

tcga_exp <- readRDS(file.path(CACHE_DIR, "tcga_exp_common.rds"))
tcga_meth_beta <- readRDS(file.path(CACHE_DIR, "tcga_meth_beta.rds"))
tcga_clin <- readRDS(file.path(CACHE_DIR, "tcga_clin_clean.rds"))
m2e <- readRDS(file.path(CACHE_DIR, "m2eqtl_discovery.rds"))

exp_features <- tcga_exp[m2e$trans_genes, , drop = FALSE]
meth_m <- beta_to_m(
  tcga_meth_beta[m2e$methylation_probes, , drop = FALSE]
)

# Author-style scaling is calculated separately for expression and methylation.
exp_scale <- scale_features_reference(exp_features)
meth_scale <- scale_features_reference(meth_m)
rownames(exp_scale$values) <- paste0("EXP__", rownames(exp_scale$values))
rownames(meth_scale$values) <- paste0("METH__", rownames(meth_scale$values))

common_samples <- Reduce(
  intersect,
  list(colnames(exp_scale$values), colnames(meth_scale$values), rownames(tcga_clin))
)
full_features <- rbind(
  exp_scale$values[, common_samples, drop = FALSE],
  meth_scale$values[, common_samples, drop = FALSE]
)

internal <- monte_carlo_m2efm(
  training_feature_by_sample = full_features,
  training_clinical = tcga_clin,
  validation_sets = list(),
  n_splits = N_MONTE_CARLO_SPLITS
)

saveRDS(internal, file.path(CACHE_DIR, "internal_full_monte_carlo.rds"))
data.table::fwrite(
  internal$metrics,
  file.path(RESULTS_DIR, "03_internal_full_monte_carlo_metrics.csv")
)

summary_metrics <- summarize_cindex_metrics(internal$metrics)
data.table::fwrite(
  summary_metrics,
  file.path(RESULTS_DIR, "03_internal_full_summary.csv")
)

save_metric_boxplot(
  internal$metrics,
  file.path(RESULTS_DIR, "03_internal_full_boxplot.pdf"),
  "TCGA internal validation: M2EFM Meth+Exp+Clin"
)

# Fit full model on all eligible matched TCGA samples.
final_internal <- fit_final_and_validate(
  training_feature_by_sample = full_features,
  training_clinical = tcga_clin
)
saveRDS(
  list(
    result = final_internal,
    expression_scaling = exp_scale[c("means", "ranges")],
    methylation_scaling = meth_scale[c("means", "ranges")],
    trans_genes = m2e$trans_genes,
    methylation_probes = m2e$methylation_probes
  ),
  file.path(RESULTS_DIR, "03_final_internal_full_model.rds")
)
data.table::fwrite(
  final_internal$metrics,
  file.path(RESULTS_DIR, "03_final_internal_full_metrics.csv")
)

message("Internal full M2EFM training completed.")
print(summary_metrics)
