#!/usr/bin/env Rscript

# TCGA OS evaluation for:
#   - M2EFM Exp
#   - M2EFM Exp+Clin
#   - rorS
#   - rorS+Clin
#   - Clin (recomputed as an audit control)
#
# Uses the same 743-patient cohort, row order, seeds, and 70/30 splits
# as results/04_internal_os_meth_exp_clin_100splits.rds.
#
# Run from repository root:
#   M2EFM_SPLITS=2   Rscript 06_run_expression_rors_os.R
#   M2EFM_SPLITS=100 Rscript 06_run_expression_rors_os.R

options(stringsAsFactors = FALSE, warn = 1)

ROOT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
R_DIR <- file.path(ROOT_DIR, "R")
DATA_DIR <- file.path(ROOT_DIR, "inst", "extdata", "csv_output")
RESULT_DIR <- file.path(ROOT_DIR, "results")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

N_SPLITS <- suppressWarnings(as.integer(Sys.getenv("M2EFM_SPLITS", "100")))
if (is.na(N_SPLITS) || N_SPLITS < 1L || N_SPLITS > 100L) {
  stop("M2EFM_SPLITS must be an integer from 1 to 100.")
}
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

summarize_result <- function(x, model) {
  x <- as.numeric(x)
  valid <- is.finite(x)
  if (!any(valid)) {
    return(data.frame(
      model = model, n_expected = length(x), n_valid = 0,
      mean = NA_real_, median = NA_real_, sd = NA_real_,
      q25 = NA_real_, q75 = NA_real_, min = NA_real_, max = NA_real_
    ))
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

same_split_lists <- function(a, b, n) {
  if (length(a) < n || length(b) < n) return(FALSE)
  all(vapply(seq_len(n), function(i) identical(as.character(a[[i]]), as.character(b[[i]])), logical(1)))
}

EXP_FILE <- file.path(DATA_DIR, "TCGA_BRCA_EXP.csv")
METH_FILE <- file.path(DATA_DIR, "TCGA_BRCA_METH.csv")
M2EQTL_FILE <- file.path(RESULT_DIR, "01_TCGA_OS_m2eQTL.rds")
REFERENCE_FILE <- file.path(RESULT_DIR, "04_internal_os_meth_exp_clin_100splits.rds")
RORS_FILE <- file.path(RESULT_DIR, "05_tcga_computed_rorS.csv")
REFERENCE_CINDEX_FILE <- file.path(RESULT_DIR, "04_internal_os_split_cindices_100splits.csv")

for (f in c(EXP_FILE, METH_FILE, M2EQTL_FILE, REFERENCE_FILE, RORS_FILE)) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}

cat("\nLoading reference 100-split evaluation...\n")
reference <- readRDS(REFERENCE_FILE)
required_reference_fields <- c("clin", "train_list", "test_list")
missing_reference_fields <- setdiff(required_reference_fields, names(reference))
if (length(missing_reference_fields)) {
  stop("Reference object lacks: ", paste(missing_reference_fields, collapse = ", "))
}

clin <- reference$clin
if (is.null(rownames(clin)) || anyDuplicated(rownames(clin))) {
  stop("Reference clinical table has invalid sample IDs.")
}
cohort_ids <- rownames(clin)

TIME_COL <- "OVERALL.SURVIVAL"
EVENT_COL <- "overall.survival.indicator"
COVARIATES <- c("pathologic_stage", "age.Dx")
needed_clin <- c(TIME_COL, EVENT_COL, COVARIATES)
missing_clin <- setdiff(needed_clin, names(clin))
if (length(missing_clin)) stop("Reference clinical table lacks: ", paste(missing_clin, collapse = ", "))

cat("Reference OS cohort: ", length(cohort_ids), " patients\n", sep = "")

m2e_os <- readRDS(M2EQTL_FILE)
expression_only_genes <- as.character(
  get_genes(m2e_os, gene_type = "trans", integrate_data = FALSE)
)
selected_probes <- as.character(get_probes(m2e_os))

cat("Expression-only genes from this run: ", length(expression_only_genes), "\n", sep = "")
cat("Selected methylation probes (loaded only for evaluate() compatibility): ",
    length(selected_probes), "\n", sep = "")

cat("\nReading TCGA expression and methylation...\n")
exp_all <- read_profile_csv(EXP_FILE, "gene_id")
meth_all <- read_profile_csv(METH_FILE, "probe_id")

exp_all <- exp_all[, grep("-01", colnames(exp_all), fixed = TRUE), drop = FALSE]
meth_all <- meth_all[, grep("-01", colnames(meth_all), fixed = TRUE), drop = FALSE]
colnames(exp_all) <- substr(colnames(exp_all), 1, 15)
colnames(meth_all) <- substr(colnames(meth_all), 1, 15)

if (anyDuplicated(colnames(exp_all))) stop("Duplicated expression IDs after truncation")
if (anyDuplicated(colnames(meth_all))) stop("Duplicated methylation IDs after truncation")

missing_genes <- setdiff(expression_only_genes, rownames(exp_all))
missing_probes <- setdiff(selected_probes, rownames(meth_all))
missing_exp_samples <- setdiff(cohort_ids, colnames(exp_all))
missing_meth_samples <- setdiff(cohort_ids, colnames(meth_all))

if (length(missing_genes)) stop("Expression-only genes missing: ", paste(missing_genes, collapse = ", "))
if (length(missing_probes)) stop("Selected probes missing: ", paste(missing_probes, collapse = ", "))
if (length(missing_exp_samples)) stop("Reference patients missing from expression: ", paste(missing_exp_samples, collapse = ", "))
if (length(missing_meth_samples)) stop("Reference patients missing from methylation: ", paste(missing_meth_samples, collapse = ", "))

exp_selected <- exp_all[expression_only_genes, cohort_ids, drop = FALSE]
meth_selected <- meth_all[selected_probes, cohort_ids, drop = FALSE]

if (any(meth_selected < 0 | meth_selected > 1)) {
  stop("Selected methylation values are not Beta-values in [0,1].")
}

exp_object <- ExpressionData$new(exp_selected)
meth_object <- MethylationData$new(meth_selected)

cat("Reading computed continuous rorS...\n")
ror <- data.table::fread(RORS_FILE, data.table = FALSE, check.names = FALSE)
if (!all(c("sample_id", "rorS") %in% names(ror))) {
  stop("rorS file must contain sample_id and rorS columns.")
}
ror$sample_id <- trimws(as.character(ror$sample_id))
ror$rorS <- suppressWarnings(as.numeric(ror$rorS))
if (anyNA(ror$sample_id) || any(!nzchar(ror$sample_id)) || anyDuplicated(ror$sample_id)) {
  stop("rorS file has invalid or duplicated sample IDs.")
}
if (anyNA(ror$rorS) || any(!is.finite(ror$rorS))) stop("rorS file has invalid scores.")

missing_rors <- setdiff(cohort_ids, ror$sample_id)
if (length(missing_rors)) stop("Reference patients missing rorS: ", paste(missing_rors, collapse = ", "))

ror <- ror[match(cohort_ids, ror$sample_id), , drop = FALSE]
pam50risk <- data.frame(rorS = ror$rorS, row.names = cohort_ids, check.names = FALSE)

pre_run_audit <- data.frame(
  metric = c(
    "reference_OS_patients",
    "OS_events",
    "OS_censored",
    "expression_only_genes",
    "rorS_scores_available",
    "rorS_min",
    "rorS_max",
    "monte_carlo_splits",
    "training_proportion",
    "ridge_alpha"
  ),
  value = c(
    nrow(clin),
    sum(clin[[EVENT_COL]] == 1),
    sum(clin[[EVENT_COL]] == 0),
    length(expression_only_genes),
    nrow(pam50risk),
    min(pam50risk$rorS),
    max(pam50risk$rorS),
    N_SPLITS,
    0.70,
    0
  )
)

cat("\nPre-run audit:\n")
print(pre_run_audit, row.names = FALSE)
write.csv(
  pre_run_audit,
  file.path(RESULT_DIR, paste0("06_expression_rors_audit_", TAG, ".csv")),
  row.names = FALSE
)

cat("\nStarting original evaluate() for expression-only and rorS models...\n")
run_start <- Sys.time()

expression_rors <- evaluate(
  meth = meth_object,
  exp = exp_object,
  clin = clin,
  eqtls = m2e_os,
  gene_type = "trans",
  covariates = COVARIATES,
  intersect_data = FALSE,
  integrate_data = FALSE,
  single_data = "exp",
  time_to_event_col = TIME_COL,
  event_col = EVENT_COL,
  rescale = TRUE,
  monte_carlo_size = N_SPLITS,
  training_prop = 0.70,
  pam50risk = pam50risk,
  alpha = 0
)

run_end <- Sys.time()
elapsed_minutes <- as.numeric(difftime(run_end, run_start, units = "mins"))

needed_result_fields <- c(
  "pred_res", "comb_res", "covar_res", "pam50_res", "pam50_comb_res",
  "train_list", "test_list"
)
missing_result_fields <- setdiff(needed_result_fields, names(expression_rors))
if (length(missing_result_fields)) {
  stop("evaluate() output lacks: ", paste(missing_result_fields, collapse = ", "))
}

same_train <- same_split_lists(expression_rors$train_list, reference$train_list, N_SPLITS)
same_test <- same_split_lists(expression_rors$test_list, reference$test_list, N_SPLITS)
if (!same_train || !same_test) {
  stop("The expression/rorS run did not reproduce the reference train/test split IDs.")
}

saveRDS(
  expression_rors,
  file.path(RESULT_DIR, paste0("06_expression_rors_os_", TAG, ".rds")),
  compress = "xz"
)

split_results <- rbind(
  data.frame(split = seq_len(N_SPLITS), model = "Clin", c_index = as.numeric(expression_rors$covar_res)),
  data.frame(split = seq_len(N_SPLITS), model = "M2EFM Exp+Clin", c_index = as.numeric(expression_rors$comb_res)),
  data.frame(split = seq_len(N_SPLITS), model = "rorS+Clin", c_index = as.numeric(expression_rors$pam50_comb_res)),
  data.frame(split = seq_len(N_SPLITS), model = "M2EFM Exp", c_index = as.numeric(expression_rors$pred_res)),
  data.frame(split = seq_len(N_SPLITS), model = "rorS", c_index = as.numeric(expression_rors$pam50_res))
)

write.csv(
  split_results,
  file.path(RESULT_DIR, paste0("06_expression_rors_split_cindices_", TAG, ".csv")),
  row.names = FALSE
)

summary_results <- rbind(
  summarize_result(expression_rors$covar_res, "Clin"),
  summarize_result(expression_rors$comb_res, "M2EFM Exp+Clin"),
  summarize_result(expression_rors$pam50_comb_res, "rorS+Clin"),
  summarize_result(expression_rors$pred_res, "M2EFM Exp"),
  summarize_result(expression_rors$pam50_res, "rorS")
)
summary_results$elapsed_minutes <- elapsed_minutes

write.csv(
  summary_results,
  file.path(RESULT_DIR, paste0("06_expression_rors_summary_", TAG, ".csv")),
  row.names = FALSE
)

# When the full 100-split run is complete, combine it with the three models
# already generated by step 04. This yields seven of the nine Figure 2A models.
if (N_SPLITS == 100L && file.exists(REFERENCE_CINDEX_FILE)) {
  old_split <- data.table::fread(REFERENCE_CINDEX_FILE, data.table = FALSE)
  new_without_duplicate_clin <- split_results[split_results$model != "Clin", , drop = FALSE]
  seven_split <- rbind(old_split, new_without_duplicate_clin)

  model_order <- c(
    "Clin",
    "M2EFM Meth+Exp+Clin",
    "M2EFM Exp+Clin",
    "rorS+Clin",
    "M2EFM Meth+Exp",
    "M2EFM Exp",
    "rorS"
  )
  seven_split$model <- factor(seven_split$model, levels = model_order)
  seven_split <- seven_split[order(seven_split$model, seven_split$split), , drop = FALSE]
  seven_split$model <- as.character(seven_split$model)

  write.csv(
    seven_split,
    file.path(RESULT_DIR, "06_seven_model_split_cindices_100splits.csv"),
    row.names = FALSE
  )

  seven_summary <- do.call(
    rbind,
    lapply(model_order, function(model_name) {
      summarize_result(
        seven_split$c_index[seven_split$model == model_name],
        model_name
      )
    })
  )

  write.csv(
    seven_summary,
    file.path(RESULT_DIR, "06_seven_model_summary_100splits.csv"),
    row.names = FALSE
  )
}

capture.output(
  {
    cat("TCGA OS expression-only and rorS evaluation\n")
    cat("Reference cohort: ", nrow(clin), " patients\n", sep = "")
    cat("Expression-only signature: ", length(expression_only_genes), " genes\n", sep = "")
    cat("Same train splits as step 04: ", same_train, "\n", sep = "")
    cat("Same test splits as step 04: ", same_test, "\n", sep = "")
    cat("Started: ", format(run_start), "\n", sep = "")
    cat("Finished: ", format(run_end), "\n", sep = "")
    cat("Elapsed minutes: ", elapsed_minutes, "\n\n", sep = "")
    print(pre_run_audit)
    print(summary_results)
    sessionInfo()
  },
  file = file.path(RESULT_DIR, paste0("06_expression_rors_session_info_", TAG, ".txt"))
)

cat("\nExpression-only and rorS summary:\n")
print(summary_results, row.names = FALSE)
cat("\nSame train splits as step 04: ", same_train, "\n", sep = "")
cat("Same test splits as step 04: ", same_test, "\n", sep = "")
cat("\nFinished successfully. Results: ", RESULT_DIR, "\n", sep = "")
