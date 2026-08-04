#!/usr/bin/env Rscript

# ============================================================
# STEP 11 — FIGURE 2B EXPRESSION-ONLY OS MODELS
#
# Input:
#   results/10B_external_os_combat_rors_ready.rds
#
# Models reconstructed (7/9 available):
#   Clin
#   M2EFM Exp+Clin
#   Cox Exp+Clin
#   rorS+Clin
#   M2EFM Exp
#   Cox Exp
#   rorS
#
# Datasets:
#   TCGA      — internal 30% test set for each split
#   Terunuma  — full external validation cohort for each split
#   Kao       — full external validation cohort for each split
#
# Unavailable:
#   NCA Exp+Clin
#   NCA Exp
# because the exact NCA gene signature is not available.
#
# Author-style workflow reproduced:
#   - full TCGA expression cohort;
#   - 100 seeds, 70/30 random splits;
#   - expression rescaled per gene using the full TCGA cohort:
#       (value - TCGA gene mean) / TCGA gene range;
#   - external cohorts transformed with those same TCGA parameters;
#   - Cox-Ridge molecular model: glmnet alpha = 0;
#   - unpenalized Cox model combines molecular risk + age + stage;
#   - clinical-only and rorS+clinical Cox models trained on each
#     TCGA training split;
#   - all models evaluated by Harrell C-index.
#
# Compatibility:
#   - cox.ties = "breslow" is set explicitly for glmnet 5.x.
#
# Usage:
#   Smoke test:
#     M2EFM_SPLITS=2 Rscript 11_run_figure2B_os_models.R
#
#   Full run:
#     M2EFM_SPLITS=100 nohup Rscript 11_run_figure2B_os_models.R \
#       > logs/11_figure2B_100splits.log 2>&1 &
#
#   Resume:
#     Re-run the same command. Completed splits are read from the
#     progress checkpoint and are not recomputed.
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1,
  width = 180
)

start_time <- Sys.time()

root_dir <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

result_dir <- file.path(
  root_dir,
  "results"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_file <- file.path(
  result_dir,
  "10B_external_os_combat_rors_ready.rds"
)

required_packages <- c(
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
    paste(
      names(package_ok)[!package_ok],
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(survcomp)
  library(matrixStats)
})

if (!file.exists(input_file)) {
  stop(
    "Missing input: ",
    input_file
  )
}

n_splits <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "M2EFM_SPLITS",
      unset = "2"
    )
  )
)

if (
  is.na(n_splits) ||
    n_splits < 1L ||
    n_splits > 100L
) {
  stop(
    "M2EFM_SPLITS must be an integer from 1 to 100."
  )
}

training_prop <- 0.70

suffix <- sprintf(
  "%03dsplits",
  n_splits
)

progress_file <- file.path(
  result_dir,
  paste0(
    "11_figure2B_progress_",
    suffix,
    ".rds"
  )
)

split_file <- file.path(
  result_dir,
  paste0(
    "11_figure2B_split_cindices_",
    suffix,
    ".csv"
  )
)

summary_file <- file.path(
  result_dir,
  paste0(
    "11_figure2B_summary_",
    suffix,
    ".csv"
  )
)

lambda_file <- file.path(
  result_dir,
  paste0(
    "11_figure2B_lambdas_",
    suffix,
    ".csv"
  )
)

sample_file <- file.path(
  result_dir,
  paste0(
    "11_figure2B_TCGA_split_sample_ids_",
    suffix,
    ".csv"
  )
)

audit_file <- file.path(
  result_dir,
  paste0(
    "11_figure2B_audit_",
    suffix,
    ".csv"
  )
)

result_file <- file.path(
  result_dir,
  paste0(
    "11_figure2B_result_",
    suffix,
    ".rds"
  )
)

session_file <- file.path(
  result_dir,
  paste0(
    "11_figure2B_sessionInfo_",
    suffix,
    ".txt"
  )
)


# ------------------------------------------------------------
# 1. Constants
# ------------------------------------------------------------

cohorts <- c(
  "TCGA",
  "Terunuma",
  "Kao"
)

model_order <- c(
  "Clin",
  "M2EFM Exp+Clin",
  "Cox Exp+Clin",
  "rorS+Clin",
  "M2EFM Exp",
  "Cox Exp",
  "rorS"
)

time_col <- "OVERALL.SURVIVAL"
event_col <- "overall.survival.indicator"
age_col <- "age.Dx"
stage_col <- "pathologic_stage"

stage_levels <- c(
  "Stage I",
  "Stage II",
  "Stage III",
  "Stage IV"
)


# ------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------

validate_expression <- function(
    x,
    cohort
) {

  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (
    is.null(rownames(x)) ||
      is.null(colnames(x))
  ) {
    stop(
      cohort,
      " expression lacks gene/sample names."
    )
  }

  if (
    anyDuplicated(rownames(x)) ||
      anyDuplicated(colnames(x))
  ) {
    stop(
      cohort,
      " expression contains duplicated genes/samples."
    )
  }

  if (
    anyNA(x) ||
      any(!is.finite(x))
  ) {
    stop(
      cohort,
      " expression contains non-finite values."
    )
  }

  x
}

prepare_clinical <- function(
    x,
    expression_sample_ids,
    cohort
) {

  if (is.null(rownames(x))) {
    stop(
      cohort,
      " clinical table lacks sample row names."
    )
  }

  required_columns <- c(
    time_col,
    event_col,
    age_col,
    stage_col,
    "rorS"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(x)
  )

  if (length(missing_columns) > 0L) {
    stop(
      cohort,
      " clinical table lacks: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  if (!identical(
    expression_sample_ids,
    rownames(x)
  )) {
    stop(
      cohort,
      " expression and clinical samples are not aligned."
    )
  }

  x[[time_col]] <- suppressWarnings(
    as.numeric(
      x[[time_col]]
    )
  )

  x[[event_col]] <- suppressWarnings(
    as.numeric(
      x[[event_col]]
    )
  )

  x[[age_col]] <- suppressWarnings(
    as.numeric(
      x[[age_col]]
    )
  )

  x$rorS <- suppressWarnings(
    as.numeric(
      x$rorS
    )
  )

  x[[stage_col]] <- factor(
    as.character(
      x[[stage_col]]
    ),
    levels = stage_levels
  )

  complete <- complete.cases(
    x[
      ,
      required_columns,
      drop = FALSE
    ]
  )

  finite <- (
    is.finite(x[[time_col]]) &
      is.finite(x[[event_col]]) &
      is.finite(x[[age_col]]) &
      is.finite(x$rorS)
  )

  valid_event <- x[[event_col]] %in% c(
    0,
    1
  )

  if (!all(
    complete &
      finite &
      valid_event
  )) {
    stop(
      cohort,
      " contains incomplete/invalid values after Step 10B."
    )
  }

  x
}

author_scale_parameters <- function(
    tcga_expression
) {

  gene_mean <- matrixStats::rowMeans2(
    tcga_expression
  )

  gene_min <- matrixStats::rowMins(
    tcga_expression
  )

  gene_max <- matrixStats::rowMaxs(
    tcga_expression
  )

  gene_range <- gene_max -
    gene_min

  names(gene_mean) <- rownames(
    tcga_expression
  )

  names(gene_range) <- rownames(
    tcga_expression
  )

  bad_range <- (
    !is.finite(gene_range) |
      gene_range <= 0
  )

  if (any(bad_range)) {
    stop(
      "Invalid TCGA gene ranges for: ",
      paste(
        head(
          names(gene_range)[bad_range],
          20L
        ),
        collapse = ", "
      )
    )
  }

  list(
    means = gene_mean,
    ranges = gene_range
  )
}

apply_author_scaling <- function(
    expression,
    means,
    ranges,
    cohort
) {

  if (!identical(
    rownames(expression),
    names(means)
  )) {
    stop(
      cohort,
      " gene order differs from scaling parameters."
    )
  }

  scaled <- sweep(
    expression,
    1L,
    means,
    FUN = "-"
  )

  scaled <- sweep(
    scaled,
    1L,
    ranges,
    FUN = "/"
  )

  if (
    anyNA(scaled) ||
      any(!is.finite(scaled))
  ) {
    stop(
      cohort,
      " scaling produced non-finite values."
    )
  }

  scaled
}

make_foldid <- function(
    n,
    nfolds = 10L
) {

  if (n < nfolds) {
    stop(
      "Training sample count is smaller than nfolds."
    )
  }

  sample(
    rep(
      seq_len(nfolds),
      length.out = n
    )
  )
}

fit_ridge_cox <- function(
    x,
    y,
    foldid
) {

  glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 0,
    foldid = foldid,
    standardize = FALSE,
    type.measure = "deviance",
    grouped = TRUE,
    cox.ties = "breslow"
  )
}

predict_ridge <- function(
    fit,
    x
) {

  prediction <- predict(
    fit,
    newx = x,
    s = "lambda.min",
    type = "response"
  )

  as.numeric(
    prediction
  )
}

fit_clinical_cox <- function(
    clinical,
    outcome
) {

  model_data <- data.frame(
    outcome = outcome,
    age.Dx = clinical[[age_col]],
    pathologic_stage = clinical[[stage_col]],
    check.names = FALSE
  )

  survival::coxph(
    outcome ~ age.Dx + pathologic_stage,
    data = model_data,
    ties = "breslow",
    x = TRUE,
    model = TRUE
  )
}

fit_combined_cox <- function(
    molecular_risk,
    clinical,
    outcome
) {

  model_data <- data.frame(
    outcome = outcome,
    pred = as.numeric(
      molecular_risk
    ),
    age.Dx = clinical[[age_col]],
    pathologic_stage = clinical[[stage_col]],
    check.names = FALSE
  )

  survival::coxph(
    outcome ~ pred + age.Dx + pathologic_stage,
    data = model_data,
    ties = "breslow",
    x = TRUE,
    model = TRUE
  )
}

fit_rors_clinical_cox <- function(
    clinical,
    outcome
) {

  model_data <- data.frame(
    outcome = outcome,
    pred = clinical$rorS,
    age.Dx = clinical[[age_col]],
    pathologic_stage = clinical[[stage_col]],
    check.names = FALSE
  )

  survival::coxph(
    outcome ~ pred + age.Dx + pathologic_stage,
    data = model_data,
    ties = "breslow",
    x = TRUE,
    model = TRUE
  )
}

clinical_newdata <- function(
    clinical
) {

  data.frame(
    age.Dx = clinical[[age_col]],
    pathologic_stage = factor(
      as.character(
        clinical[[stage_col]]
      ),
      levels = stage_levels
    ),
    check.names = FALSE
  )
}

combined_newdata <- function(
    molecular_risk,
    clinical
) {

  data.frame(
    pred = as.numeric(
      molecular_risk
    ),
    age.Dx = clinical[[age_col]],
    pathologic_stage = factor(
      as.character(
        clinical[[stage_col]]
      ),
      levels = stage_levels
    ),
    check.names = FALSE
  )
}

rors_newdata <- function(
    clinical
) {

  data.frame(
    pred = clinical$rorS,
    age.Dx = clinical[[age_col]],
    pathologic_stage = factor(
      as.character(
        clinical[[stage_col]]
      ),
      levels = stage_levels
    ),
    check.names = FALSE
  )
}

predict_cox_lp <- function(
    fit,
    newdata
) {

  as.numeric(
    predict(
      fit,
      newdata = newdata,
      type = "lp",
      reference = "sample"
    )
  )
}

c_index <- function(
    risk,
    clinical
) {

  result <- survcomp::concordance.index(
    x = as.numeric(risk),
    surv.time = clinical[[time_col]],
    surv.event = clinical[[event_col]],
    method = "noether"
  )

  value <- as.numeric(
    result$c.index
  )

  if (
    length(value) != 1L ||
      is.na(value) ||
      !is.finite(value)
  ) {
    stop(
      "Invalid C-index returned."
    )
  }

  value
}

append_cindex_rows <- function(
    split,
    dataset,
    values
) {

  if (!identical(
    names(values),
    model_order
  )) {
    stop(
      "Model result order mismatch."
    )
  }

  data.frame(
    split = split,
    dataset = dataset,
    model = names(values),
    c_index = as.numeric(values),
    stringsAsFactors = FALSE
  )
}

summarize_cindices <- function(
    x
) {

  groups <- split(
    x,
    interaction(
      x$dataset,
      factor(
        x$model,
        levels = model_order
      ),
      drop = TRUE,
      lex.order = TRUE
    )
  )

  summary_rows <- lapply(
    groups,
    function(group) {

      values <- group$c_index

      data.frame(
        dataset = group$dataset[[1L]],
        model = group$model[[1L]],
        n_expected = n_splits,
        n_valid = sum(
          is.finite(values)
        ),
        mean = mean(
          values,
          na.rm = TRUE
        ),
        median = median(
          values,
          na.rm = TRUE
        ),
        sd = stats::sd(
          values,
          na.rm = TRUE
        ),
        q25 = unname(
          quantile(
            values,
            0.25,
            na.rm = TRUE
          )
        ),
        q75 = unname(
          quantile(
            values,
            0.75,
            na.rm = TRUE
          )
        ),
        min = min(
          values,
          na.rm = TRUE
        ),
        max = max(
          values,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  output <- do.call(
    rbind,
    summary_rows
  )

  output$dataset <- factor(
    output$dataset,
    levels = cohorts
  )

  output$model <- factor(
    output$model,
    levels = model_order
  )

  output <- output[
    order(
      output$dataset,
      output$model
    ),
    ,
    drop = FALSE
  ]

  output$dataset <- as.character(
    output$dataset
  )

  output$model <- as.character(
    output$model
  )

  rownames(output) <- NULL

  output
}


# ------------------------------------------------------------
# 3. Load and validate model-ready data
# ------------------------------------------------------------

cat(
  "Loading ",
  input_file,
  "...\n",
  sep = ""
)

ready <- readRDS(
  input_file
)

missing_cohorts <- setdiff(
  cohorts,
  names(ready)
)

if (length(missing_cohorts) > 0L) {
  stop(
    "Input lacks cohorts: ",
    paste(
      missing_cohorts,
      collapse = ", "
    )
  )
}

expression <- lapply(
  cohorts,
  function(cohort) {
    validate_expression(
      ready[[cohort]]$expression,
      cohort
    )
  }
)

names(expression) <- cohorts

clinical <- lapply(
  cohorts,
  function(cohort) {
    prepare_clinical(
      ready[[cohort]]$clinical,
      colnames(
        expression[[cohort]]
      ),
      cohort
    )
  }
)

names(clinical) <- cohorts

common_genes <- rownames(
  expression$TCGA
)

for (cohort in cohorts) {
  if (!identical(
    rownames(
      expression[[cohort]]
    ),
    common_genes
  )) {
    stop(
      "Common gene order differs for ",
      cohort,
      "."
    )
  }
}

signature_genes <- as.character(
  ready$external_expression_signature
)

signature_genes <- unique(
  signature_genes
)

if (length(signature_genes) != 113L) {
  stop(
    "Expected 113 M2EFM genes but found ",
    length(signature_genes),
    "."
  )
}

missing_signature <- setdiff(
  signature_genes,
  common_genes
)

if (length(missing_signature) > 0L) {
  stop(
    "M2EFM genes missing from common expression matrix: ",
    paste(
      missing_signature,
      collapse = ", "
    )
  )
}

input_md5 <- unname(
  tools::md5sum(
    input_file
  )
)


# ------------------------------------------------------------
# 4. Apply the repository's full-TCGA mean/range scaling
# ------------------------------------------------------------

cat(
  "Applying author-style TCGA mean/range scaling...\n"
)

scale_parameters <- author_scale_parameters(
  expression$TCGA
)

expression_scaled <- lapply(
  cohorts,
  function(cohort) {
    apply_author_scaling(
      expression[[cohort]],
      scale_parameters$means,
      scale_parameters$ranges,
      cohort
    )
  }
)

names(expression_scaled) <- cohorts

tcga_m2efm <- t(
  expression_scaled$TCGA[
    signature_genes,
    ,
    drop = FALSE
  ]
)

terunuma_m2efm <- t(
  expression_scaled$Terunuma[
    signature_genes,
    ,
    drop = FALSE
  ]
)

kao_m2efm <- t(
  expression_scaled$Kao[
    signature_genes,
    ,
    drop = FALSE
  ]
)

tcga_cox <- t(
  expression_scaled$TCGA
)

terunuma_cox <- t(
  expression_scaled$Terunuma
)

kao_cox <- t(
  expression_scaled$Kao
)

stopifnot(
  identical(
    rownames(tcga_m2efm),
    rownames(clinical$TCGA)
  ),
  identical(
    rownames(terunuma_m2efm),
    rownames(clinical$Terunuma)
  ),
  identical(
    rownames(kao_m2efm),
    rownames(clinical$Kao)
  ),
  identical(
    rownames(tcga_cox),
    rownames(clinical$TCGA)
  ),
  identical(
    rownames(terunuma_cox),
    rownames(clinical$Terunuma)
  ),
  identical(
    rownames(kao_cox),
    rownames(clinical$Kao)
  )
)


# ------------------------------------------------------------
# 5. Initialize or resume progress
# ------------------------------------------------------------

progress <- list(
  input_md5 = input_md5,
  n_splits = n_splits,
  training_prop = training_prop,
  model_order = model_order,
  completed_splits = integer(0),
  cindices = data.frame(),
  lambdas = data.frame(),
  split_samples = data.frame(),
  started_at = start_time,
  updated_at = start_time
)

if (file.exists(progress_file)) {

  saved_progress <- readRDS(
    progress_file
  )

  compatible <- (
    identical(
      saved_progress$input_md5,
      input_md5
    ) &&
      identical(
        saved_progress$n_splits,
        n_splits
      ) &&
      identical(
        saved_progress$training_prop,
        training_prop
      ) &&
      identical(
        saved_progress$model_order,
        model_order
      )
  )

  if (!compatible) {
    stop(
      "Existing progress checkpoint is incompatible with current data/settings: ",
      progress_file
    )
  }

  progress <- saved_progress

  cat(
    "Resuming from ",
    length(
      progress$completed_splits
    ),
    " completed split(s).\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 6. Run Monte-Carlo splits
# ------------------------------------------------------------

n_tcga <- nrow(
  tcga_cox
)

train_size <- as.integer(
  training_prop *
    n_tcga
)

all_indices <- seq_len(
  n_tcga
)

external_data <- list(
  Terunuma = list(
    m2efm = terunuma_m2efm,
    cox = terunuma_cox,
    clinical = clinical$Terunuma
  ),
  Kao = list(
    m2efm = kao_m2efm,
    cox = kao_cox,
    clinical = clinical$Kao
  )
)

for (split_number in seq_len(n_splits)) {

  if (split_number %in% progress$completed_splits) {
    next
  }

  split_start <- Sys.time()

  cat(
    "\nSeed ",
    split_number,
    "\n",
    sep = ""
  )

  set.seed(
    split_number
  )

  train_index <- sample(
    all_indices,
    size = train_size,
    replace = FALSE
  )

  test_index <- setdiff(
    all_indices,
    train_index
  )

  train_ids <- rownames(
    tcga_cox
  )[train_index]

  test_ids <- rownames(
    tcga_cox
  )[test_index]

  clin_train <- clinical$TCGA[
    train_ids,
    ,
    drop = FALSE
  ]

  clin_test <- clinical$TCGA[
    test_ids,
    ,
    drop = FALSE
  ]

  outcome_train <- survival::Surv(
    clin_train[[time_col]],
    clin_train[[event_col]]
  )

  foldid <- make_foldid(
    length(train_ids),
    nfolds = 10L
  )

  # M2EFM expression-only Cox-Ridge.
  m2efm_fit <- fit_ridge_cox(
    x = tcga_m2efm[
      train_ids,
      ,
      drop = FALSE
    ],
    y = outcome_train,
    foldid = foldid
  )

  m2efm_train_risk <- predict_ridge(
    m2efm_fit,
    tcga_m2efm[
      train_ids,
      ,
      drop = FALSE
    ]
  )

  m2efm_final_fit <- fit_combined_cox(
    molecular_risk = m2efm_train_risk,
    clinical = clin_train,
    outcome = outcome_train
  )

  # Transcriptome-wide Cox-Ridge comparison.
  cox_fit <- fit_ridge_cox(
    x = tcga_cox[
      train_ids,
      ,
      drop = FALSE
    ],
    y = outcome_train,
    foldid = foldid
  )

  cox_train_risk <- predict_ridge(
    cox_fit,
    tcga_cox[
      train_ids,
      ,
      drop = FALSE
    ]
  )

  cox_final_fit <- fit_combined_cox(
    molecular_risk = cox_train_risk,
    clinical = clin_train,
    outcome = outcome_train
  )

  # Clinical and PAM50 rorS comparison models.
  clinical_fit <- fit_clinical_cox(
    clinical = clin_train,
    outcome = outcome_train
  )

  rors_clinical_fit <- fit_rors_clinical_cox(
    clinical = clin_train,
    outcome = outcome_train
  )

  # TCGA internal predictions.
  tcga_m2efm_risk <- predict_ridge(
    m2efm_fit,
    tcga_m2efm[
      test_ids,
      ,
      drop = FALSE
    ]
  )

  tcga_cox_risk <- predict_ridge(
    cox_fit,
    tcga_cox[
      test_ids,
      ,
      drop = FALSE
    ]
  )

  tcga_values <- c(
    "Clin" = c_index(
      predict_cox_lp(
        clinical_fit,
        clinical_newdata(
          clin_test
        )
      ),
      clin_test
    ),
    "M2EFM Exp+Clin" = c_index(
      predict_cox_lp(
        m2efm_final_fit,
        combined_newdata(
          tcga_m2efm_risk,
          clin_test
        )
      ),
      clin_test
    ),
    "Cox Exp+Clin" = c_index(
      predict_cox_lp(
        cox_final_fit,
        combined_newdata(
          tcga_cox_risk,
          clin_test
        )
      ),
      clin_test
    ),
    "rorS+Clin" = c_index(
      predict_cox_lp(
        rors_clinical_fit,
        rors_newdata(
          clin_test
        )
      ),
      clin_test
    ),
    "M2EFM Exp" = c_index(
      tcga_m2efm_risk,
      clin_test
    ),
    "Cox Exp" = c_index(
      tcga_cox_risk,
      clin_test
    ),
    "rorS" = c_index(
      clin_test$rorS,
      clin_test
    )
  )

  split_rows <- list(
    append_cindex_rows(
      split_number,
      "TCGA",
      tcga_values
    )
  )

  # Full external validation predictions.
  for (external_cohort in c(
    "Terunuma",
    "Kao"
  )) {

    current <- external_data[
      [external_cohort]
    ]

    ext_m2efm_risk <- predict_ridge(
      m2efm_fit,
      current$m2efm
    )

    ext_cox_risk <- predict_ridge(
      cox_fit,
      current$cox
    )

    ext_clin <- current$clinical

    external_values <- c(
      "Clin" = c_index(
        predict_cox_lp(
          clinical_fit,
          clinical_newdata(
            ext_clin
          )
        ),
        ext_clin
      ),
      "M2EFM Exp+Clin" = c_index(
        predict_cox_lp(
          m2efm_final_fit,
          combined_newdata(
            ext_m2efm_risk,
            ext_clin
          )
        ),
        ext_clin
      ),
      "Cox Exp+Clin" = c_index(
        predict_cox_lp(
          cox_final_fit,
          combined_newdata(
            ext_cox_risk,
            ext_clin
          )
        ),
        ext_clin
      ),
      "rorS+Clin" = c_index(
        predict_cox_lp(
          rors_clinical_fit,
          rors_newdata(
            ext_clin
          )
        ),
        ext_clin
      ),
      "M2EFM Exp" = c_index(
        ext_m2efm_risk,
        ext_clin
      ),
      "Cox Exp" = c_index(
        ext_cox_risk,
        ext_clin
      ),
      "rorS" = c_index(
        ext_clin$rorS,
        ext_clin
      )
    )

    split_rows[
      [length(split_rows) + 1L]
    ] <- append_cindex_rows(
      split_number,
      external_cohort,
      external_values
    )
  }

  progress$cindices <- rbind(
    progress$cindices,
    do.call(
      rbind,
      split_rows
    )
  )

  progress$lambdas <- rbind(
    progress$lambdas,
    data.frame(
      split = split_number,
      M2EFM_lambda_min = m2efm_fit$lambda.min,
      M2EFM_lambda_1se = m2efm_fit$lambda.1se,
      Cox_lambda_min = cox_fit$lambda.min,
      Cox_lambda_1se = cox_fit$lambda.1se,
      stringsAsFactors = FALSE
    )
  )

  progress$split_samples <- rbind(
    progress$split_samples,
    data.frame(
      split = split_number,
      set = c(
        rep(
          "train",
          length(train_ids)
        ),
        rep(
          "test",
          length(test_ids)
        )
      ),
      sample_id = c(
        train_ids,
        test_ids
      ),
      stringsAsFactors = FALSE
    )
  )

  progress$completed_splits <- sort(
    unique(
      c(
        progress$completed_splits,
        split_number
      )
    )
  )

  progress$updated_at <- Sys.time()

  saveRDS(
    progress,
    progress_file,
    compress = "xz"
  )

  write.csv(
    progress$cindices,
    split_file,
    row.names = FALSE
  )

  write.csv(
    progress$lambdas,
    lambda_file,
    row.names = FALSE
  )

  write.csv(
    progress$split_samples,
    sample_file,
    row.names = FALSE
  )

  split_minutes <- as.numeric(
    difftime(
      Sys.time(),
      split_start,
      units = "mins"
    )
  )

  current_tcga_combined <- tcga_values[
    "M2EFM Exp+Clin"
  ]

  running_mean <- mean(
    progress$cindices$c_index[
      progress$cindices$dataset ==
        "TCGA" &
        progress$cindices$model ==
          "M2EFM Exp+Clin"
    ]
  )

  cat(
    "TCGA: ",
    paste(
      names(tcga_values),
      sprintf(
        "%.4f",
        tcga_values
      ),
      sep = "=",
      collapse = " | "
    ),
    "\n",
    sep = ""
  )

  cat(
    "Split minutes: ",
    round(
      split_minutes,
      3
    ),
    " | Running mean M2EFM Exp+Clin: ",
    round(
      running_mean,
      6
    ),
    "\n",
    sep = ""
  )

  rm(
    m2efm_fit,
    cox_fit,
    m2efm_final_fit,
    cox_final_fit,
    clinical_fit,
    rors_clinical_fit
  )

  invisible(
    gc()
  )
}


# ------------------------------------------------------------
# 7. Validate and summarize
# ------------------------------------------------------------

expected_rows <- (
  n_splits *
    length(cohorts) *
    length(model_order)
)

if (nrow(progress$cindices) != expected_rows) {
  stop(
    "Expected ",
    expected_rows,
    " C-index rows but found ",
    nrow(progress$cindices),
    "."
  )
}

combination_counts <- table(
  progress$cindices$dataset,
  progress$cindices$model
)

if (any(combination_counts != n_splits)) {
  stop(
    "At least one dataset/model combination lacks ",
    n_splits,
    " results."
  )
}

summary_table <- summarize_cindices(
  progress$cindices
)

write.csv(
  summary_table,
  summary_file,
  row.names = FALSE
)

audit <- data.frame(
  metric = c(
    "input_file",
    "input_md5",
    "requested_splits",
    "completed_splits",
    "training_proportion",
    "TCGA_samples",
    "TCGA_train_samples_per_split",
    "TCGA_test_samples_per_split",
    "TCGA_OS_events",
    "Terunuma_samples",
    "Terunuma_OS_events",
    "Kao_samples",
    "Kao_OS_events",
    "common_expression_genes",
    "M2EFM_external_signature_genes",
    "models_reconstructed",
    "expected_C_index_values",
    "observed_C_index_values",
    "glmnet_alpha",
    "glmnet_standardize",
    "glmnet_cox_ties",
    "expression_scaling",
    "NCA_models_available"
  ),
  value = c(
    normalizePath(
      input_file,
      winslash = "/",
      mustWork = TRUE
    ),
    input_md5,
    n_splits,
    length(
      progress$completed_splits
    ),
    training_prop,
    nrow(
      clinical$TCGA
    ),
    train_size,
    nrow(
      clinical$TCGA
    ) -
      train_size,
    sum(
      clinical$TCGA[[event_col]] ==
        1
    ),
    nrow(
      clinical$Terunuma
    ),
    sum(
      clinical$Terunuma[[event_col]] ==
        1
    ),
    nrow(
      clinical$Kao
    ),
    sum(
      clinical$Kao[[event_col]] ==
        1
    ),
    length(
      common_genes
    ),
    length(
      signature_genes
    ),
    length(
      model_order
    ),
    expected_rows,
    nrow(
      progress$cindices
    ),
    0,
    FALSE,
    "breslow",
    "(value - full-TCGA gene mean) / full-TCGA gene range",
    FALSE
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  audit_file,
  row.names = FALSE
)

end_time <- Sys.time()

elapsed_minutes <- as.numeric(
  difftime(
    end_time,
    start_time,
    units = "mins"
  )
)

final_result <- list(
  settings = list(
    n_splits = n_splits,
    training_prop = training_prop,
    seeds = seq_len(
      n_splits
    ),
    models = model_order,
    datasets = cohorts,
    time_column = time_col,
    event_column = event_col,
    age_column = age_col,
    stage_column = stage_col,
    stage_levels = stage_levels,
    glmnet_alpha = 0,
    glmnet_standardize = FALSE,
    glmnet_cox_ties = "breslow",
    expression_scaling = paste(
      "(value - full-TCGA gene mean)",
      "/ full-TCGA gene range"
    ),
    NCA_available = FALSE
  ),
  cindices = progress$cindices,
  summary = summary_table,
  lambdas = progress$lambdas,
  split_samples = progress$split_samples,
  audit = audit,
  input_md5 = input_md5,
  elapsed_minutes_this_invocation = elapsed_minutes,
  finished_at = end_time
)

saveRDS(
  final_result,
  result_file,
  compress = "xz"
)

capture.output(
  {
    cat("STEP 11 — FIGURE 2B EXPRESSION-ONLY OS MODELS\n")
    cat("Finished: ", format(end_time), "\n", sep = "")
    cat("Elapsed minutes this invocation: ", elapsed_minutes, "\n\n", sep = "")

    cat("SUMMARY\n")
    print(summary_table)

    cat("\nAUDIT\n")
    print(audit)

    cat("\nSESSION INFO\n")
    sessionInfo()
  },
  file = session_file
)


# ------------------------------------------------------------
# 8. Console report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 11 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nFIGURE 2B SUMMARY — 7 AVAILABLE MODELS\n")
print(
  summary_table[
    ,
    c(
      "dataset",
      "model",
      "n_valid",
      "mean",
      "median",
      "sd"
    )
  ],
  row.names = FALSE
)

cat(
  "\nC-index values saved:\n",
  split_file,
  "\n",
  sep = ""
)

cat(
  "\nSummary saved:\n",
  summary_file,
  "\n",
  sep = ""
)

cat(
  "\nResult object saved:\n",
  result_file,
  "\n",
  sep = ""
)

cat(
  "\nElapsed minutes this invocation: ",
  round(
    elapsed_minutes,
    3
  ),
  "\n",
  sep = ""
)

cat("\nFinished successfully.\n")
