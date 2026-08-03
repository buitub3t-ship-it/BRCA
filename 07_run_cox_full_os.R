#!/usr/bin/env Rscript

# TCGA OS Cox-Ridge baselines for the two remaining Figure 2A models:
#   - Cox Meth+Exp
#   - Cox Meth+Exp+Clin
#
# This script uses this replication's available data:
#   - all available expression genes
#   - the most variable available methylation probes by MAD, up to the
#     number of available expression genes
#
# In the current processed files this is expected to be:
#   10,000 expression genes + 7,313 available methylation probes.
#
# The original paper used no gene/probe outcome pre-selection for this
# baseline, but MAD-filtered methylation probes to roughly the number of
# genes. Because this replication's methylation CSV contains only 7,313
# processed probes, all 7,313 are retained when fewer probes than genes
# are available.
#
# Uses the same 743-patient cohort, seeds, and 70/30 splits as step 04.
#
# Run from repository root:
#   cd ~/BRCA
#   M2EFM_SPLITS=2 Rscript 07_run_cox_full_os.R
#   M2EFM_SPLITS=100 Rscript 07_run_cox_full_os.R

options(stringsAsFactors = FALSE, warn = 1)

ROOT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
R_DIR <- file.path(ROOT_DIR, "R")
DATA_DIR <- file.path(ROOT_DIR, "inst", "extdata", "csv_output")
RESULT_DIR <- file.path(ROOT_DIR, "results")

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

N_SPLITS <- suppressWarnings(as.integer(Sys.getenv("M2EFM_SPLITS", "100")))
if (is.na(N_SPLITS) || N_SPLITS < 1L) {
  stop("M2EFM_SPLITS must be a positive integer.")
}

TAG <- sprintf("%03dsplits", N_SPLITS)

cat("Repository root: ", ROOT_DIR, "\n", sep = "")
cat("Monte-Carlo splits: ", N_SPLITS, "\n", sep = "")

required_packages <- c(
  "R6",
  "data.table",
  "impute",
  "lumi",
  "glmnet",
  "survival",
  "survcomp",
  "matrixStats"
)

package_ok <- vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

if (!all(package_ok)) {
  stop(
    "Missing packages: ",
    paste(names(package_ok)[!package_ok], collapse = ", ")
  )
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

# Lock the first-stage Cox glmnet fit to the glmnet 5.0 behavior used
# in steps 04 and 06, and silence the announced future default change.
cv.glmnet_original <- glmnet::cv.glmnet
cv.glmnet <- function(...) {
  cv.glmnet_original(..., cox.ties = "breslow")
}

for (f in c(
  "ProfileData.R",
  "MethylationData.R",
  "ExpressionData.R",
  "m2efm.R"
)) {
  path <- file.path(R_DIR, f)
  if (!file.exists(path)) {
    stop("Missing source file: ", path)
  }
  source(path, local = .GlobalEnv)
}

read_profile_csv <- function(file, id_column) {
  if (!file.exists(file)) {
    stop("Missing input file: ", file)
  }

  x <- data.table::fread(
    file,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c("NA", "NaN", "")
  )

  if (!id_column %in% names(x)) {
    stop(basename(file), " lacks column ", id_column)
  }

  ids <- trimws(as.character(x[[id_column]]))

  if (anyNA(ids) || any(!nzchar(ids))) {
    stop(basename(file), " has empty feature IDs")
  }

  if (anyDuplicated(ids)) {
    stop(basename(file), " has duplicated feature IDs")
  }

  x[[id_column]] <- NULL
  rownames(x) <- ids

  x[] <- lapply(
    x,
    function(v) suppressWarnings(as.numeric(v))
  )

  if (anyNA(x)) {
    stop(
      basename(file),
      " contains NA or non-numeric molecular values"
    )
  }

  if (anyDuplicated(colnames(x))) {
    stop(basename(file), " has duplicated sample IDs")
  }

  x
}

summarize_result <- function(x, model) {
  x <- as.numeric(x)
  valid <- is.finite(x)

  if (!any(valid)) {
    return(
      data.frame(
        model = model,
        n_expected = length(x),
        n_valid = 0,
        mean = NA_real_,
        median = NA_real_,
        sd = NA_real_,
        q25 = NA_real_,
        q75 = NA_real_,
        min = NA_real_,
        max = NA_real_
      )
    )
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
  if (length(a) < n || length(b) < n) {
    return(FALSE)
  }

  all(
    vapply(
      seq_len(n),
      function(i) {
        identical(
          as.character(a[[i]]),
          as.character(b[[i]])
        )
      },
      logical(1)
    )
  )
}

EXP_FILE <- file.path(DATA_DIR, "TCGA_BRCA_EXP.csv")
METH_FILE <- file.path(DATA_DIR, "TCGA_BRCA_METH.csv")

REFERENCE_FILE <- file.path(
  RESULT_DIR,
  "04_internal_os_meth_exp_clin_100splits.rds"
)

SEVEN_MODEL_FILE <- file.path(
  RESULT_DIR,
  "06_seven_model_split_cindices_100splits.csv"
)

for (f in c(EXP_FILE, METH_FILE, REFERENCE_FILE)) {
  if (!file.exists(f)) {
    stop("Missing input file: ", f)
  }
}

cat("\nLoading reference 100-split evaluation...\n")

reference <- readRDS(REFERENCE_FILE)

required_reference_fields <- c(
  "clin",
  "covar_res",
  "train_list",
  "test_list"
)

missing_reference_fields <- setdiff(
  required_reference_fields,
  names(reference)
)

if (length(missing_reference_fields)) {
  stop(
    "Reference object lacks: ",
    paste(missing_reference_fields, collapse = ", ")
  )
}

clin <- reference$clin

if (
  is.null(rownames(clin)) ||
  anyDuplicated(rownames(clin))
) {
  stop("Reference clinical table has invalid sample IDs.")
}

cohort_ids <- rownames(clin)

TIME_COL <- "OVERALL.SURVIVAL"
EVENT_COL <- "overall.survival.indicator"
COVARIATES <- c("pathologic_stage", "age.Dx")

needed_clin <- c(
  TIME_COL,
  EVENT_COL,
  COVARIATES
)

missing_clin <- setdiff(
  needed_clin,
  names(clin)
)

if (length(missing_clin)) {
  stop(
    "Reference clinical table lacks: ",
    paste(missing_clin, collapse = ", ")
  )
}

cat(
  "Reference OS cohort: ",
  length(cohort_ids),
  " patients\n",
  sep = ""
)

cat("\nReading complete available expression and methylation profiles...\n")

exp_all <- read_profile_csv(
  EXP_FILE,
  "gene_id"
)

meth_all <- read_profile_csv(
  METH_FILE,
  "probe_id"
)

exp_all <- exp_all[
  ,
  grep("-01", colnames(exp_all), fixed = TRUE),
  drop = FALSE
]

meth_all <- meth_all[
  ,
  grep("-01", colnames(meth_all), fixed = TRUE),
  drop = FALSE
]

colnames(exp_all) <- substr(
  colnames(exp_all),
  1,
  15
)

colnames(meth_all) <- substr(
  colnames(meth_all),
  1,
  15
)

if (anyDuplicated(colnames(exp_all))) {
  stop("Duplicated expression IDs after truncation")
}

if (anyDuplicated(colnames(meth_all))) {
  stop("Duplicated methylation IDs after truncation")
}

missing_exp_samples <- setdiff(
  cohort_ids,
  colnames(exp_all)
)

missing_meth_samples <- setdiff(
  cohort_ids,
  colnames(meth_all)
)

if (length(missing_exp_samples)) {
  stop(
    "Reference patients missing from expression: ",
    paste(missing_exp_samples, collapse = ", ")
  )
}

if (length(missing_meth_samples)) {
  stop(
    "Reference patients missing from methylation: ",
    paste(missing_meth_samples, collapse = ", ")
  )
}

exp_full <- exp_all[
  ,
  cohort_ids,
  drop = FALSE
]

meth_beta_full <- meth_all[
  ,
  cohort_ids,
  drop = FALSE
]

if (
  any(meth_beta_full < 0) ||
  any(meth_beta_full > 1)
) {
  stop("Methylation input is not Beta-value data in [0,1].")
}

cat(
  "Available expression: ",
  nrow(exp_full),
  " genes x ",
  ncol(exp_full),
  " patients\n",
  sep = ""
)

cat(
  "Available methylation: ",
  nrow(meth_beta_full),
  " probes x ",
  ncol(meth_beta_full),
  " patients\n",
  sep = ""
)

cat("\nCalculating methylation MAD on M-values...\n")

meth_for_mad <- MethylationData$new(
  meth_beta_full
)

meth_m_full <- meth_for_mad$m_values

if (any(!is.finite(as.matrix(meth_m_full)))) {
  stop("Non-finite methylation M-values detected.")
}

meth_mad <- matrixStats::rowMads(
  as.matrix(meth_m_full),
  na.rm = TRUE
)

names(meth_mad) <- rownames(meth_m_full)

if (any(!is.finite(meth_mad))) {
  stop("Non-finite methylation MAD values detected.")
}

n_expression_genes <- nrow(exp_full)

n_target_probes <- min(
  n_expression_genes,
  nrow(meth_beta_full)
)

selected_probes <- names(
  sort(
    meth_mad,
    decreasing = TRUE
  )
)[seq_len(n_target_probes)]

if (length(selected_probes) != n_target_probes) {
  stop("Failed to select the requested number of probes.")
}

cat(
  "Cox baseline features: ",
  n_expression_genes,
  " expression genes + ",
  length(selected_probes),
  " methylation probes = ",
  n_expression_genes + length(selected_probes),
  " molecular predictors\n",
  sep = ""
)

if (nrow(meth_beta_full) < n_expression_genes) {
  cat(
    "Note: fewer processed methylation probes than expression genes are ",
    "available, so all available probes are retained.\n",
    sep = ""
  )
}

signature_file <- file.path(
  RESULT_DIR,
  paste0(
    "07_cox_available_feature_signature_",
    TAG,
    ".txt"
  )
)

writeLines(
  c(
    rownames(exp_full),
    selected_probes
  ),
  signature_file
)

probe_mad_table <- data.frame(
  probe_id = selected_probes,
  MAD_M_value = as.numeric(
    meth_mad[selected_probes]
  ),
  rank = seq_along(selected_probes),
  stringsAsFactors = FALSE
)

write.csv(
  probe_mad_table,
  file.path(
    RESULT_DIR,
    paste0(
      "07_cox_selected_probe_mad_",
      TAG,
      ".csv"
    )
  ),
  row.names = FALSE
)

exp_object <- ExpressionData$new(
  exp_full
)

meth_object <- MethylationData$new(
  meth_beta_full[selected_probes, , drop = FALSE]
)

pre_run_audit <- data.frame(
  metric = c(
    "reference_OS_patients",
    "OS_events",
    "OS_censored",
    "available_expression_genes",
    "available_methylation_probes",
    "selected_methylation_probes_by_MAD",
    "total_molecular_predictors",
    "methylation_probe_target_rule",
    "monte_carlo_splits",
    "training_proportion",
    "ridge_alpha",
    "glmnet_cox_ties"
  ),
  value = c(
    nrow(clin),
    sum(clin[[EVENT_COL]] == 1),
    sum(clin[[EVENT_COL]] == 0),
    n_expression_genes,
    nrow(meth_beta_full),
    length(selected_probes),
    n_expression_genes + length(selected_probes),
    paste0(
      "min(expression genes, available probes) = ",
      n_target_probes
    ),
    N_SPLITS,
    0.70,
    0,
    "breslow"
  ),
  stringsAsFactors = FALSE
)

cat("\nPre-run audit:\n")
print(pre_run_audit, row.names = FALSE)

write.csv(
  pre_run_audit,
  file.path(
    RESULT_DIR,
    paste0(
      "07_cox_full_audit_",
      TAG,
      ".csv"
    )
  ),
  row.names = FALSE
)

cat("\nStarting Cox-Ridge Meth+Exp evaluation...\n")
cat("This run is expected to be substantially slower than steps 04 and 06.\n")

run_start <- Sys.time()

cox_full <- evaluate(
  meth = meth_object,
  exp = exp_object,
  clin = clin,
  gene_type = "pre",
  gene_sig = signature_file,
  num_genes = n_expression_genes,
  num_probes = length(selected_probes),
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

elapsed_minutes <- as.numeric(
  difftime(
    run_end,
    run_start,
    units = "mins"
  )
)

required_result_fields <- c(
  "pred_res",
  "comb_res",
  "covar_res",
  "train_list",
  "test_list"
)

missing_result_fields <- setdiff(
  required_result_fields,
  names(cox_full)
)

if (length(missing_result_fields)) {
  stop(
    "evaluate() output lacks: ",
    paste(missing_result_fields, collapse = ", ")
  )
}

same_train <- same_split_lists(
  cox_full$train_list,
  reference$train_list,
  N_SPLITS
)

same_test <- same_split_lists(
  cox_full$test_list,
  reference$test_list,
  N_SPLITS
)

if (!same_train || !same_test) {
  stop(
    "The Cox run did not reproduce the reference train/test split IDs."
  )
}

clinical_difference <- max(
  abs(
    as.numeric(cox_full$covar_res) -
      as.numeric(reference$covar_res[seq_len(N_SPLITS)])
  ),
  na.rm = TRUE
)

if (
  !is.finite(clinical_difference) ||
  clinical_difference > 1e-10
) {
  stop(
    "Clinical-only control differs from step 04. Maximum difference: ",
    clinical_difference
  )
}

saveRDS(
  cox_full,
  file.path(
    RESULT_DIR,
    paste0(
      "07_cox_full_os_",
      TAG,
      ".rds"
    )
  ),
  compress = "xz"
)

split_results <- rbind(
  data.frame(
    split = seq_len(N_SPLITS),
    model = "Cox Meth+Exp+Clin",
    c_index = as.numeric(cox_full$comb_res)
  ),
  data.frame(
    split = seq_len(N_SPLITS),
    model = "Cox Meth+Exp",
    c_index = as.numeric(cox_full$pred_res)
  )
)

write.csv(
  split_results,
  file.path(
    RESULT_DIR,
    paste0(
      "07_cox_full_split_cindices_",
      TAG,
      ".csv"
    )
  ),
  row.names = FALSE
)

summary_results <- rbind(
  summarize_result(
    cox_full$comb_res,
    "Cox Meth+Exp+Clin"
  ),
  summarize_result(
    cox_full$pred_res,
    "Cox Meth+Exp"
  )
)

summary_results$elapsed_minutes <- elapsed_minutes

write.csv(
  summary_results,
  file.path(
    RESULT_DIR,
    paste0(
      "07_cox_full_summary_",
      TAG,
      ".csv"
    )
  ),
  row.names = FALSE
)

if (
  N_SPLITS == 100L &&
  file.exists(SEVEN_MODEL_FILE)
) {
  seven_split <- data.table::fread(
    SEVEN_MODEL_FILE,
    data.table = FALSE
  )

  nine_split <- rbind(
    seven_split,
    split_results
  )

  model_order <- c(
    "Clin",
    "M2EFM Meth+Exp+Clin",
    "M2EFM Exp+Clin",
    "Cox Meth+Exp+Clin",
    "rorS+Clin",
    "M2EFM Meth+Exp",
    "M2EFM Exp",
    "Cox Meth+Exp",
    "rorS"
  )

  nine_split$model <- factor(
    nine_split$model,
    levels = model_order
  )

  nine_split <- nine_split[
    order(
      nine_split$model,
      nine_split$split
    ),
    ,
    drop = FALSE
  ]

  nine_split$model <- as.character(
    nine_split$model
  )

  write.csv(
    nine_split,
    file.path(
      RESULT_DIR,
      "07_nine_model_split_cindices_100splits.csv"
    ),
    row.names = FALSE
  )

  nine_summary <- do.call(
    rbind,
    lapply(
      model_order,
      function(model_name) {
        summarize_result(
          nine_split$c_index[
            nine_split$model == model_name
          ],
          model_name
        )
      }
    )
  )

  write.csv(
    nine_summary,
    file.path(
      RESULT_DIR,
      "07_nine_model_summary_100splits.csv"
    ),
    row.names = FALSE
  )
}

capture.output(
  {
    cat("TCGA OS available-data Cox-Ridge baseline\n")
    cat(
      "Expression genes: ",
      n_expression_genes,
      "\n",
      sep = ""
    )
    cat(
      "Selected methylation probes: ",
      length(selected_probes),
      "\n",
      sep = ""
    )
    cat(
      "Total molecular predictors: ",
      n_expression_genes + length(selected_probes),
      "\n",
      sep = ""
    )
    cat(
      "Same train splits as step 04: ",
      same_train,
      "\n",
      sep = ""
    )
    cat(
      "Same test splits as step 04: ",
      same_test,
      "\n",
      sep = ""
    )
    cat(
      "Maximum clinical control difference: ",
      clinical_difference,
      "\n",
      sep = ""
    )
    cat(
      "Started: ",
      format(run_start),
      "\n",
      sep = ""
    )
    cat(
      "Finished: ",
      format(run_end),
      "\n",
      sep = ""
    )
    cat(
      "Elapsed minutes: ",
      elapsed_minutes,
      "\n\n",
      sep = ""
    )

    print(pre_run_audit)
    print(summary_results)
    sessionInfo()
  },
  file = file.path(
    RESULT_DIR,
    paste0(
      "07_cox_full_session_info_",
      TAG,
      ".txt"
    )
  )
)

cat("\nCox-Ridge summary:\n")
print(summary_results, row.names = FALSE)

cat(
  "\nSame train splits as step 04: ",
  same_train,
  "\n",
  sep = ""
)

cat(
  "Same test splits as step 04: ",
  same_test,
  "\n",
  sep = ""
)

cat(
  "Maximum clinical control difference: ",
  clinical_difference,
  "\n",
  sep = ""
)

cat(
  "\nFinished successfully. Results: ",
  RESULT_DIR,
  "\n",
  sep = ""
)
