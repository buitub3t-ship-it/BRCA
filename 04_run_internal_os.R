#!/usr/bin/env Rscript

# Internal TCGA OS evaluation using this replication's own M2EFM features.
# Run from the repository root, for example:
#   cd ~/BRCA
#   M2EFM_SPLITS=2 Rscript 04_run_internal_os.R

options(stringsAsFactors = FALSE, warn = 1)

ROOT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
R_DIR <- file.path(ROOT_DIR, "R")
DATA_DIR <- file.path(ROOT_DIR, "inst", "extdata", "csv_output")
RESULT_DIR <- file.path(ROOT_DIR, "results")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

N_SPLITS <- suppressWarnings(as.integer(Sys.getenv("M2EFM_SPLITS", "100")))
if (is.na(N_SPLITS) || N_SPLITS < 1L) stop("M2EFM_SPLITS must be a positive integer.")
TAG <- sprintf("%03dsplits", N_SPLITS)

cat("Repository root: ", ROOT_DIR, "\n", sep = "")
cat("Monte-Carlo splits: ", N_SPLITS, "\n", sep = "")

required_packages <- c(
  "R6", "data.table", "impute", "lumi", "glmnet",
  "survival", "survcomp", "matrixStats"
)
package_ok <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(package_ok)) {
  stop("Missing packages: ", paste(names(package_ok)[!package_ok], collapse = ", "))
}

suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(impute)
  library(lumi)
  library(glmnet)
  library(survival)
  library(survcomp)
  library(matrixStats)
})

# Original author classes and model functions.
for (f in c("ProfileData.R", "MethylationData.R", "ExpressionData.R", "m2efm.R")) {
  path <- file.path(R_DIR, f)
  if (!file.exists(path)) stop("Missing source file: ", path)
  source(path, local = .GlobalEnv)
}

read_profile_csv <- function(file, id_column) {
  if (!file.exists(file)) stop("Missing input file: ", file)
  x <- data.table::fread(
    file,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c("NA", "NaN", "")
  )
  if (!id_column %in% names(x)) stop(basename(file), " lacks column ", id_column)
  ids <- trimws(as.character(x[[id_column]]))
  if (anyNA(ids) || any(!nzchar(ids))) stop(basename(file), " has empty feature IDs")
  if (anyDuplicated(ids)) stop(basename(file), " has duplicated feature IDs")
  x[[id_column]] <- NULL
  rownames(x) <- ids
  x[] <- lapply(x, function(v) suppressWarnings(as.numeric(v)))
  if (anyNA(x)) stop(basename(file), " contains NA or non-numeric molecular values")
  if (anyDuplicated(colnames(x))) stop(basename(file), " has duplicated sample IDs")
  x
}

EXP_FILE <- file.path(DATA_DIR, "TCGA_BRCA_EXP.csv")
METH_FILE <- file.path(DATA_DIR, "TCGA_BRCA_METH.csv")
CLIN_FILE <- file.path(DATA_DIR, "TCGA_BRCA_CLIN.csv")
M2EQTL_FILE <- file.path(RESULT_DIR, "01_TCGA_OS_m2eQTL.rds")

for (f in c(EXP_FILE, METH_FILE, CLIN_FILE, M2EQTL_FILE)) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}

m2e_os <- readRDS(M2EQTL_FILE)
trans_genes <- as.character(get_genes(m2e_os, "trans", integrate_data = TRUE))
selected_probes <- as.character(get_probes(m2e_os))
expression_only_genes <- as.character(
  get_genes(m2e_os, "trans", integrate_data = FALSE)
)

feature_manifest <- data.frame(
  metric = c(
    "trans_m2eGenes", "selected_methylation_probes",
    "integrated_features", "expression_only_genes"
  ),
  value = c(
    length(trans_genes), length(selected_probes),
    length(trans_genes) + length(selected_probes),
    length(expression_only_genes)
  )
)
cat("\nFeature manifest from this run:\n")
print(feature_manifest, row.names = FALSE)
write.csv(
  feature_manifest,
  file.path(RESULT_DIR, paste0("04_internal_os_feature_manifest_", TAG, ".csv")),
  row.names = FALSE
)

cat("\nReading molecular data...\n")
tcga_exp <- read_profile_csv(EXP_FILE, "gene_id")
tcga_meth <- read_profile_csv(METH_FILE, "probe_id")

cat("Reading clinical data...\n")
tcga_clin <- data.table::fread(
  CLIN_FILE,
  data.table = FALSE,
  check.names = FALSE,
  na.strings = c("NA", "N/A", "--", "")
)
if (!"sample_id" %in% names(tcga_clin)) stop("Clinical file lacks sample_id")
clinical_ids <- trimws(as.character(tcga_clin$sample_id))
if (anyNA(clinical_ids) || any(!nzchar(clinical_ids))) stop("Empty clinical sample IDs")
if (anyDuplicated(clinical_ids)) stop("Duplicated clinical sample IDs")
rownames(tcga_clin) <- clinical_ids
tcga_clin$sample_id <- NULL

TIME_COL <- "OVERALL.SURVIVAL"
EVENT_COL <- "overall.survival.indicator"
COVARIATES <- c("pathologic_stage", "age.Dx")
needed_clin <- c(TIME_COL, EVENT_COL, COVARIATES)
missing_clin <- setdiff(needed_clin, names(tcga_clin))
if (length(missing_clin)) stop("Missing clinical columns: ", paste(missing_clin, collapse = ", "))

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
tcga_clin[[TIME_COL]] <- as_num(tcga_clin[[TIME_COL]])
tcga_clin[[EVENT_COL]] <- as_num(tcga_clin[[EVENT_COL]])
tcga_clin$age.Dx <- as_num(tcga_clin$age.Dx)
tcga_clin$pathologic_stage <- factor(
  reduce_stage(as.character(tcga_clin$pathologic_stage)),
  levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
)

# Primary tumors only; preserve the author's 15-character TCGA ID convention.
tcga_exp <- tcga_exp[, grep("-01", colnames(tcga_exp), fixed = TRUE), drop = FALSE]
tcga_meth <- tcga_meth[, grep("-01", colnames(tcga_meth), fixed = TRUE), drop = FALSE]
colnames(tcga_exp) <- substr(colnames(tcga_exp), 1, 15)
colnames(tcga_meth) <- substr(colnames(tcga_meth), 1, 15)
if (anyDuplicated(colnames(tcga_exp))) stop("Duplicated expression IDs after truncation")
if (anyDuplicated(colnames(tcga_meth))) stop("Duplicated methylation IDs after truncation")

common_samples <- Reduce(
  intersect,
  list(colnames(tcga_exp), colnames(tcga_meth), rownames(tcga_clin))
)
if (!length(common_samples)) stop("No common TCGA samples across expression, methylation, and clinical data")

clin_common <- tcga_clin[common_samples, , drop = FALSE]
complete <- complete.cases(clin_common[, needed_clin, drop = FALSE]) &
  is.finite(clin_common[[TIME_COL]]) &
  is.finite(clin_common[[EVENT_COL]]) &
  is.finite(clin_common$age.Dx)
usable_samples <- rownames(clin_common)[complete]
if (length(usable_samples) < 20L) stop("Too few complete OS samples: ", length(usable_samples))

tcga_clin <- clin_common[usable_samples, , drop = FALSE]
tcga_clin$pathologic_stage <- droplevels(tcga_clin$pathologic_stage)
if (!all(unique(tcga_clin[[EVENT_COL]]) %in% c(0, 1))) stop("OS event is not coded 0/1")
if (any(tcga_clin[[TIME_COL]] < 0)) stop("Negative OS time detected")

missing_genes <- setdiff(trans_genes, rownames(tcga_exp))
missing_probes <- setdiff(selected_probes, rownames(tcga_meth))
if (length(missing_genes)) stop("Selected expression genes missing: ", paste(missing_genes, collapse = ", "))
if (length(missing_probes)) stop("Selected methylation probes missing: ", paste(missing_probes, collapse = ", "))

exp_selected <- tcga_exp[trans_genes, usable_samples, drop = FALSE]
meth_selected <- tcga_meth[selected_probes, usable_samples, drop = FALSE]
if (any(meth_selected < 0 | meth_selected > 1)) stop("Methylation input is not Beta-value data in [0,1]")

exp_object <- ExpressionData$new(exp_selected)
meth_object <- MethylationData$new(meth_selected)

pre_run_audit <- data.frame(
  metric = c(
    "common_samples_before_complete_case_filter",
    "samples_removed_missing_OS_or_covariates",
    "usable_OS_samples", "OS_events", "OS_censored",
    "trans_m2eGenes", "selected_methylation_probes",
    "integrated_features", "monte_carlo_splits",
    "training_proportion", "ridge_alpha"
  ),
  value = c(
    length(common_samples), sum(!complete), length(usable_samples),
    sum(tcga_clin[[EVENT_COL]] == 1), sum(tcga_clin[[EVENT_COL]] == 0),
    length(trans_genes), length(selected_probes),
    length(trans_genes) + length(selected_probes),
    N_SPLITS, 0.70, 0
  )
)
cat("\nPre-run audit:\n")
print(pre_run_audit, row.names = FALSE)
write.csv(
  pre_run_audit,
  file.path(RESULT_DIR, paste0("04_internal_os_pre_run_audit_", TAG, ".csv")),
  row.names = FALSE
)

cat("\nStarting original evaluate()...\n")
run_start <- Sys.time()
internal_os <- evaluate(
  meth = meth_object,
  exp = exp_object,
  clin = tcga_clin,
  eqtls = m2e_os,
  gene_type = "trans",
  covariates = COVARIATES,
  intersect_data = TRUE,
  integrate_data = TRUE,
  time_to_event_col = TIME_COL,
  event_col = EVENT_COL,
  rescale = TRUE,
  monte_carlo_size = N_SPLITS,
  training_prop = 0.70,
  alpha = 0
)
run_end <- Sys.time()
elapsed_minutes <- as.numeric(difftime(run_end, run_start, units = "mins"))

saveRDS(
  internal_os,
  file.path(RESULT_DIR, paste0("04_internal_os_meth_exp_clin_", TAG, ".rds")),
  compress = "xz"
)

split_results <- rbind(
  data.frame(split = seq_along(internal_os$covar_res), model = "Clin", c_index = as.numeric(internal_os$covar_res)),
  data.frame(split = seq_along(internal_os$comb_res), model = "M2EFM Meth+Exp+Clin", c_index = as.numeric(internal_os$comb_res)),
  data.frame(split = seq_along(internal_os$pred_res), model = "M2EFM Meth+Exp", c_index = as.numeric(internal_os$pred_res))
)
write.csv(
  split_results,
  file.path(RESULT_DIR, paste0("04_internal_os_split_cindices_", TAG, ".csv")),
  row.names = FALSE
)

summarize_result <- function(x, model) {
  x <- as.numeric(x)
  valid <- is.finite(x)
  if (!any(valid)) {
    return(data.frame(model = model, n_expected = length(x), n_valid = 0,
                      mean = NA, median = NA, sd = NA, q25 = NA, q75 = NA,
                      min = NA, max = NA))
  }
  y <- x[valid]
  data.frame(
    model = model,
    n_expected = length(x),
    n_valid = length(y),
    mean = mean(y),
    median = median(y),
    sd = if (length(y) > 1L) stats::sd(y) else NA_real_,
    q25 = unname(stats::quantile(y, 0.25)),
    q75 = unname(stats::quantile(y, 0.75)),
    min = min(y),
    max = max(y)
  )
}

summary_results <- rbind(
  summarize_result(internal_os$covar_res, "Clin"),
  summarize_result(internal_os$comb_res, "M2EFM Meth+Exp+Clin"),
  summarize_result(internal_os$pred_res, "M2EFM Meth+Exp")
)
summary_results$elapsed_minutes <- elapsed_minutes
write.csv(
  summary_results,
  file.path(RESULT_DIR, paste0("04_internal_os_summary_", TAG, ".csv")),
  row.names = FALSE
)

capture.output(
  {
    cat("Internal TCGA OS evaluation\n")
    cat("Features: ", length(trans_genes), " genes + ", length(selected_probes),
        " probes = ", length(trans_genes) + length(selected_probes), "\n", sep = "")
    cat("Started: ", format(run_start), "\n", sep = "")
    cat("Finished: ", format(run_end), "\n", sep = "")
    cat("Elapsed minutes: ", elapsed_minutes, "\n\n", sep = "")
    print(pre_run_audit)
    print(summary_results)
    sessionInfo()
  },
  file = file.path(RESULT_DIR, paste0("04_internal_os_session_info_", TAG, ".txt"))
)

cat("\nInternal OS summary:\n")
print(summary_results, row.names = FALSE)
cat("\nFinished successfully. Results: ", RESULT_DIR, "\n", sep = "")
