#!/usr/bin/env Rscript

# =====================================================================
# STEPS 14–17 — ALL REMAINING OVERALL-SURVIVAL OUTPUTS
#
# Runs in one command:
#   Step 14: Figure S1-style 5-year calibration + IPCW Brier/IBS
#   Step 15: Kao ER / ER-negative / triple-negative subgroup Table 2
#   Step 16: GO-BP pathway heatmap resembling Figure 6
#   Step 17: Author-style OS tables and output manifest
#
# IMPORTANT LIMITATIONS
# ---------------------
# 1. Table 2:
#    - Uses Kao receptor columns when they exist.
#    - Otherwise infers ER/PR/HER2 from ESR1/PGR/ERBB2 expression
#      using a two-component Gaussian mixture (or k-means fallback).
#    - The paper's exact receptor probe set was not publicly specified,
#      so inferred results are a methodological approximation.
#
# 2. Figure 6:
#    - Repeats GO Biological Process enrichment with the installed
#      Bioconductor annotation database.
#    - It cannot exactly reproduce WebGestalt 2017 database contents.
#
# 3. Table S5/S6:
#    - Copies/standardizes results already created in prior steps.
#    - NCA models and 1,000 random-gene bootstraps remain unavailable
#      unless separately supplied and run.
#
# Main input:
#   results/10B_external_os_combat_rors_ready.rds
#
# Run:
#   time Rscript 14_17_run_all_remaining_OS.R \
#     2>&1 | tee logs/14_17_run_all_remaining_OS.log
# =====================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1,
  width = 200
)

start_time <- Sys.time()

root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

results_dir <- file.path(
  root,
  "results"
)

output_dir <- file.path(
  results_dir,
  "14_17_OS_remaining"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_file <- file.path(
  results_dir,
  "10B_external_os_combat_rors_ready.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Missing input file: ",
    input_file
  )
}

core_packages <- c(
  "glmnet",
  "survival",
  "matrixStats",
  "survcomp"
)

core_ok <- vapply(
  core_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

if (!all(core_ok)) {
  stop(
    "Missing required packages: ",
    paste(
      names(core_ok)[!core_ok],
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(matrixStats)
  library(survcomp)
})

cohorts <- c(
  "TCGA",
  "Terunuma",
  "Kao"
)

stage_levels <- c(
  "Stage I",
  "Stage II",
  "Stage III",
  "Stage IV"
)

risk_levels <- c(
  "Low",
  "Medium",
  "High"
)

time_col <- "OVERALL.SURVIVAL"
event_col <- "overall.survival.indicator"
age_col <- "age.Dx"
stage_col <- "pathologic_stage"

days_per_month <- 365.25 / 12
five_year_days <- 5 * 365.25
ten_year_days <- 10 * 365.25
seed_value <- 1L

time_multipliers <- c(
  TCGA = 1,
  Terunuma = days_per_month,
  Kao = 365.25
)

status_rows <- list()
status_index <- 1L

record_status <- function(
    step,
    component,
    status,
    note = ""
) {

  status_rows[[status_index]] <<- data.frame(
    step = step,
    component = component,
    status = status,
    note = note,
    stringsAsFactors = FALSE
  )

  status_index <<- status_index + 1L

  invisible(
    NULL
  )
}


# =====================================================================
# COMMON HELPERS
# =====================================================================

validate_expression <- function(
    x,
    cohort
) {

  x <- as.matrix(
    x
  )

  storage.mode(
    x
  ) <- "double"

  if (
    is.null(
      rownames(x)
    ) ||
      is.null(
        colnames(x)
      )
  ) {
    stop(
      cohort,
      " expression lacks gene/sample names."
    )
  }

  if (
    anyDuplicated(
      rownames(x)
    ) ||
      anyDuplicated(
        colnames(x)
      )
  ) {
    stop(
      cohort,
      " expression has duplicated names."
    )
  }

  if (
    anyNA(x) ||
      any(
        !is.finite(x)
      )
  ) {
    stop(
      cohort,
      " expression has non-finite values."
    )
  }

  x
}


prepare_clinical <- function(
    x,
    sample_ids,
    cohort
) {

  required_columns <- c(
    time_col,
    event_col,
    age_col,
    stage_col
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(x)
  )

  if (length(
    missing_columns
  ) > 0L) {
    stop(
      cohort,
      " clinical table lacks: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  if (
    is.null(
      rownames(x)
    ) ||
      !identical(
        rownames(x),
        sample_ids
      )
  ) {
    stop(
      cohort,
      " clinical/expression samples are not aligned."
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

  x[[stage_col]] <- factor(
    as.character(
      x[[stage_col]]
    ),
    levels = stage_levels
  )

  valid <- complete.cases(
    x[
      ,
      required_columns,
      drop = FALSE
    ]
  ) &
    is.finite(
      x[[time_col]]
    ) &
    is.finite(
      x[[event_col]]
    ) &
    is.finite(
      x[[age_col]]
    ) &
    x[[time_col]] >= 0 &
    x[[event_col]] %in% c(
      0,
      1
    )

  if (!all(
    valid
  )) {
    stop(
      cohort,
      " contains incomplete/invalid OS values."
    )
  }

  x[[time_col]] <-
    x[[time_col]] *
    time_multipliers[[cohort]]

  x
}


scale_parameters <- function(
    tcga_expression
) {

  means <- matrixStats::rowMeans2(
    tcga_expression
  )

  ranges <- matrixStats::rowMaxs(
    tcga_expression
  ) -
    matrixStats::rowMins(
      tcga_expression
    )

  names(means) <- rownames(
    tcga_expression
  )

  names(ranges) <- rownames(
    tcga_expression
  )

  if (any(
    !is.finite(ranges) |
      ranges <= 0
  )) {
    stop(
      "Invalid TCGA expression ranges."
    )
  }

  list(
    means = means,
    ranges = ranges
  )
}


apply_scaling <- function(
    x,
    parameters,
    cohort
) {

  if (!identical(
    rownames(x),
    names(
      parameters$means
    )
  )) {
    stop(
      cohort,
      " gene order differs from scaling parameters."
    )
  }

  z <- sweep(
    x,
    1L,
    parameters$means,
    FUN = "-"
  )

  z <- sweep(
    z,
    1L,
    parameters$ranges,
    FUN = "/"
  )

  if (
    anyNA(z) ||
      any(
        !is.finite(z)
      )
  ) {
    stop(
      cohort,
      " scaling produced non-finite values."
    )
  }

  z
}


make_foldid <- function(
    n,
    folds = 10L
) {

  sample(
    rep(
      seq_len(
        folds
      ),
      length.out = n
    )
  )
}


predict_ridge_response <- function(
    fit,
    x
) {

  as.numeric(
    predict(
      fit,
      newx = x,
      s = "lambda.min",
      type = "response"
    )
  )
}


make_integrated_data <- function(
    molecular_risk,
    clinical,
    include_outcome = FALSE
) {

  output <- data.frame(
    molecular_risk = as.numeric(
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

  if (include_outcome) {
    output$survival_outcome <- survival::Surv(
      clinical[[time_col]],
      clinical[[event_col]]
    )
  }

  output
}


predict_integrated_lp <- function(
    fit,
    molecular_risk,
    clinical
) {

  newdata <- make_integrated_data(
    molecular_risk,
    clinical,
    include_outcome = FALSE
  )

  as.numeric(
    predict(
      fit,
      newdata = newdata,
      type = "lp",
      reference = "sample"
    )
  )
}


predict_event_probability <- function(
    fit,
    molecular_risk,
    clinical,
    times
) {

  newdata <- make_integrated_data(
    molecular_risk,
    clinical,
    include_outcome = FALSE
  )

  lp_zero <- as.numeric(
    predict(
      fit,
      newdata = newdata,
      type = "lp",
      reference = "zero"
    )
  )

  baseline <- survival::basehaz(
    fit,
    centered = FALSE
  )

  cumulative_hazard <- vapply(
    times,
    function(current_time) {

      index <- max(
        which(
          baseline$time <=
            current_time
        ),
        0L
      )

      if (
        length(index) == 0L ||
          index == 0L
      ) {
        return(
          0
        )
      }

      baseline$hazard[[index]]
    },
    numeric(1)
  )

  probability_matrix <- vapply(
    cumulative_hazard,
    function(h0) {
      1 -
        exp(
          -h0 *
            exp(
              lp_zero
            )
        )
    },
    numeric(
      length(
        lp_zero
      )
    )
  )

  if (is.null(
    dim(
      probability_matrix
    )
  )) {
    probability_matrix <- matrix(
      probability_matrix,
      ncol = length(
        times
      )
    )
  }

  rownames(
    probability_matrix
  ) <- rownames(
    clinical
  )

  colnames(
    probability_matrix
  ) <- as.character(
    times
  )

  probability_matrix
}


assign_risk_group <- function(
    score
) {

  cutoffs <- unname(
    quantile(
      score,
      probs = c(
        0.25,
        0.75
      ),
      type = 7
    )
  )

  group <- ifelse(
    score < cutoffs[[1L]],
    "Low",
    ifelse(
      score > cutoffs[[2L]],
      "High",
      "Medium"
    )
  )

  factor(
    group,
    levels = risk_levels
  )
}


cindex_table_row <- function(
    score,
    clinical,
    cohort,
    subgroup
) {

  if (
    length(score) < 5L ||
      sum(
        clinical[[event_col]] ==
          1
      ) < 2L
  ) {
    return(
      data.frame(
        cohort = cohort,
        subgroup = subgroup,
        samples = length(
          score
        ),
        events = sum(
          clinical[[event_col]] ==
            1
        ),
        c_index = NA_real_,
        lower_95 = NA_real_,
        upper_95 = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  result <- survcomp::concordance.index(
    x = as.numeric(
      score
    ),
    surv.time = clinical[[time_col]],
    surv.event = clinical[[event_col]],
    method = "noether"
  )

  data.frame(
    cohort = cohort,
    subgroup = subgroup,
    samples = length(
      score
    ),
    events = sum(
      clinical[[event_col]] ==
        1
    ),
    c_index = as.numeric(
      result$c.index
    ),
    lower_95 = as.numeric(
      result$lower
    ),
    upper_95 = as.numeric(
      result$upper
    ),
    stringsAsFactors = FALSE
  )
}


safe_quantile_groups <- function(
    x,
    groups = 4L
) {

  ranks <- rank(
    x,
    ties.method = "first"
  )

  factor(
    pmin(
      groups,
      ceiling(
        ranks /
          (
            length(
              ranks
            ) /
              groups
          )
      )
    ),
    levels = seq_len(
      groups
    )
  )
}


observed_event_probability <- function(
    time,
    event,
    target_time
) {

  fit <- survival::survfit(
    survival::Surv(
      time,
      event
    ) ~ 1
  )

  result <- summary(
    fit,
    times = target_time,
    extend = TRUE
  )

  1 -
    as.numeric(
      result$surv[[1L]]
    )
}


km_step_survival <- function(
    fit,
    times,
    left_limit = FALSE
) {

  adjusted_times <- times

  if (left_limit) {
    adjusted_times <-
      adjusted_times -
      pmax(
        abs(
          adjusted_times
        ),
        1
      ) *
      1e-10
  }

  indices <- findInterval(
    adjusted_times,
    fit$time
  )

  survival_values <- c(
    1,
    fit$surv
  )

  survival_values[
    indices +
      1L
  ]
}


ipcw_brier <- function(
    time,
    event,
    predicted_event_probability,
    evaluation_time
) {

  censoring_fit <- survival::survfit(
    survival::Surv(
      time,
      1 -
        event
    ) ~ 1
  )

  g_at_t <- km_step_survival(
    censoring_fit,
    evaluation_time,
    left_limit = FALSE
  )

  g_at_time_minus <- km_step_survival(
    censoring_fit,
    time,
    left_limit = TRUE
  )

  epsilon <- 1e-6

  g_at_t <- pmax(
    g_at_t,
    epsilon
  )

  g_at_time_minus <- pmax(
    g_at_time_minus,
    epsilon
  )

  event_before_t <- (
    time <=
      evaluation_time &
      event ==
        1
  )

  event_free_at_t <- (
    time >
      evaluation_time
  )

  loss <- numeric(
    length(
      time
    )
  )

  loss[event_before_t] <-
    (
      1 -
        predicted_event_probability[event_before_t]
    )^2 /
    g_at_time_minus[event_before_t]

  loss[event_free_at_t] <-
    predicted_event_probability[event_free_at_t]^2 /
    g_at_t

  mean(
    loss
  )
}


trapezoid_integral <- function(
    x,
    y
) {

  if (length(
    x
  ) < 2L) {
    return(
      NA_real_
    )
  }

  sum(
    diff(
      x
    ) *
      (
        head(
          y,
          -1L
        ) +
          tail(
            y,
            -1L
          )
      ) /
      2
  )
}


normalize_binary_status <- function(
    values
) {

  text <- tolower(
    trimws(
      as.character(
        values
      )
    )
  )

  output <- rep(
    NA_character_,
    length(
      text
    )
  )

  positive <- grepl(
    "positive|pos|yes|true|present|amplified|^1$",
    text
  )

  negative <- grepl(
    "negative|neg|no|false|absent|not amplified|^0$",
    text
  )

  output[positive] <- "Positive"
  output[negative] <- "Negative"

  output
}


find_first_column <- function(
    column_names,
    exact_candidates,
    regex_candidates
) {

  lower_names <- tolower(
    column_names
  )

  for (candidate in exact_candidates) {

    index <- which(
      lower_names ==
        tolower(
          candidate
        )
    )

    if (length(
      index
    ) > 0L) {
      return(
        column_names[
          index[[1L]]
        ]
      )
    }
  }

  for (pattern in regex_candidates) {

    index <- grep(
      pattern,
      lower_names,
      perl = TRUE
    )

    if (length(
      index
    ) > 0L) {
      return(
        column_names[
          index[[1L]]
        ]
      )
    }
  }

  NA_character_
}


infer_two_component_status <- function(
    expression_values,
    sample_ids,
    label
) {

  expression_values <- as.numeric(
    expression_values
  )

  if (
    anyNA(
      expression_values
    ) ||
      any(
        !is.finite(
          expression_values
        )
      )
  ) {
    stop(
      label,
      " expression contains invalid values."
    )
  }

  method <- "kmeans_2_components"

  if (requireNamespace(
    "mclust",
    quietly = TRUE
  )) {

    fitted <- try(
      mclust::Mclust(
        expression_values,
        G = 2,
        verbose = FALSE
      ),
      silent = TRUE
    )

    if (!inherits(
      fitted,
      "try-error"
    )) {

      cluster <- fitted$classification
      method <- "mclust_G2"

    } else {

      fitted <- stats::kmeans(
        expression_values,
        centers = 2,
        nstart = 100
      )

      cluster <- fitted$cluster
    }

  } else {

    fitted <- stats::kmeans(
      expression_values,
      centers = 2,
      nstart = 100
    )

    cluster <- fitted$cluster
  }

  cluster_means <- tapply(
    expression_values,
    cluster,
    mean
  )

  positive_cluster <- as.integer(
    names(
      which.max(
        cluster_means
      )
    )
  )

  status <- ifelse(
    cluster ==
      positive_cluster,
    "Positive",
    "Negative"
  )

  data.frame(
    sample_id = sample_ids,
    marker = label,
    expression = expression_values,
    cluster = cluster,
    status = status,
    inference_method = method,
    stringsAsFactors = FALSE
  )
}


format_cox_table <- function(
    fit
) {

  fit_summary <- summary(
    fit
  )

  data.frame(
    term = rownames(
      fit_summary$coefficients
    ),
    coefficient = fit_summary$coefficients[
      ,
      "coef"
    ],
    hazard_ratio = fit_summary$conf.int[
      ,
      "exp(coef)"
    ],
    lower_95 = fit_summary$conf.int[
      ,
      "lower .95"
    ],
    upper_95 = fit_summary$conf.int[
      ,
      "upper .95"
    ],
    p_value = fit_summary$coefficients[
      ,
      "Pr(>|z|)"
    ],
    stringsAsFactors = FALSE
  )
}


# =====================================================================
# LOAD AND PREPARE DATA
# =====================================================================

cat(
  "Loading model-ready OS data...\n"
)

ready <- readRDS(
  input_file
)

missing_cohorts <- setdiff(
  cohorts,
  names(
    ready
  )
)

if (length(
  missing_cohorts
) > 0L) {
  stop(
    "Input lacks cohorts: ",
    paste(
      missing_cohorts,
      collapse = ", "
    )
  )
}

expression <- setNames(
  lapply(
    cohorts,
    function(cohort) {
      validate_expression(
        ready[[cohort]]$expression,
        cohort
      )
    }
  ),
  cohorts
)

clinical <- setNames(
  lapply(
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
  ),
  cohorts
)

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
      cohort,
      " common-gene order differs."
    )
  }
}

signature_genes <- unique(
  as.character(
    ready$external_expression_signature
  )
)

if (length(
  signature_genes
) != 113L) {
  warning(
    "Expected 113 external-common signature genes; found ",
    length(
      signature_genes
    ),
    "."
  )
}

missing_signature <- setdiff(
  signature_genes,
  common_genes
)

if (length(
  missing_signature
) > 0L) {
  stop(
    "Signature genes missing from common expression matrix: ",
    paste(
      missing_signature,
      collapse = ", "
    )
  )
}

parameters <- scale_parameters(
  expression$TCGA
)

scaled_expression <- setNames(
  lapply(
    cohorts,
    function(cohort) {
      apply_scaling(
        expression[[cohort]],
        parameters,
        cohort
      )
    }
  ),
  cohorts
)

m2efm_x <- setNames(
  lapply(
    cohorts,
    function(cohort) {
      t(
        scaled_expression[[cohort]][
          signature_genes,
          ,
          drop = FALSE
        ]
      )
    }
  ),
  cohorts
)

cox_x <- setNames(
  lapply(
    cohorts,
    function(cohort) {
      t(
        scaled_expression[[cohort]]
      )
    }
  ),
  cohorts
)


# =====================================================================
# FIT OR LOAD FULL-TCGA M2EFM AND COX MODELS
# =====================================================================

model_checkpoint <- file.path(
  output_dir,
  "14_17_full_TCGA_models.rds"
)

input_md5 <- unname(
  tools::md5sum(
    input_file
  )
)

models <- NULL

if (file.exists(
  model_checkpoint
)) {

  checkpoint <- try(
    readRDS(
      model_checkpoint
    ),
    silent = TRUE
  )

  if (
    !inherits(
      checkpoint,
      "try-error"
    ) &&
      identical(
        checkpoint$input_md5,
        input_md5
      ) &&
      identical(
        checkpoint$signature_genes,
        signature_genes
      )
  ) {

    models <- checkpoint

    cat(
      "Loaded compatible full-model checkpoint.\n"
    )
  }
}

if (is.null(
  models
)) {

  cat(
    "Fitting full-TCGA M2EFM and transcriptome-wide Cox-Ridge models...\n"
  )

  set.seed(
    seed_value
  )

  tcga_outcome <- survival::Surv(
    clinical$TCGA[[time_col]],
    clinical$TCGA[[event_col]]
  )

  foldid <- make_foldid(
    nrow(
      m2efm_x$TCGA
    ),
    folds = 10L
  )

  m2efm_fit <- glmnet::cv.glmnet(
    x = m2efm_x$TCGA,
    y = tcga_outcome,
    family = "cox",
    alpha = 0,
    foldid = foldid,
    standardize = FALSE,
    type.measure = "deviance",
    grouped = TRUE,
    cox.ties = "breslow"
  )

  cox_fit <- glmnet::cv.glmnet(
    x = cox_x$TCGA,
    y = tcga_outcome,
    family = "cox",
    alpha = 0,
    foldid = foldid,
    standardize = FALSE,
    type.measure = "deviance",
    grouped = TRUE,
    cox.ties = "breslow"
  )

  m2efm_molecular <- setNames(
    lapply(
      cohorts,
      function(cohort) {
        predict_ridge_response(
          m2efm_fit,
          m2efm_x[[cohort]]
        )
      }
    ),
    cohorts
  )

  cox_molecular <- setNames(
    lapply(
      cohorts,
      function(cohort) {
        predict_ridge_response(
          cox_fit,
          cox_x[[cohort]]
        )
      }
    ),
    cohorts
  )

  m2efm_integrated_fit <- survival::coxph(
    survival_outcome ~
      molecular_risk +
      age.Dx +
      pathologic_stage,
    data = make_integrated_data(
      m2efm_molecular$TCGA,
      clinical$TCGA,
      include_outcome = TRUE
    ),
    ties = "breslow",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )

  cox_integrated_fit <- survival::coxph(
    survival_outcome ~
      molecular_risk +
      age.Dx +
      pathologic_stage,
    data = make_integrated_data(
      cox_molecular$TCGA,
      clinical$TCGA,
      include_outcome = TRUE
    ),
    ties = "breslow",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )

  m2efm_final_risk <- setNames(
    lapply(
      cohorts,
      function(cohort) {
        predict_integrated_lp(
          m2efm_integrated_fit,
          m2efm_molecular[[cohort]],
          clinical[[cohort]]
        )
      }
    ),
    cohorts
  )

  cox_final_risk <- setNames(
    lapply(
      cohorts,
      function(cohort) {
        predict_integrated_lp(
          cox_integrated_fit,
          cox_molecular[[cohort]],
          clinical[[cohort]]
        )
      }
    ),
    cohorts
  )

  models <- list(
    input_md5 = input_md5,
    signature_genes = signature_genes,
    parameters = parameters,
    foldid = foldid,
    m2efm_fit = m2efm_fit,
    cox_fit = cox_fit,
    m2efm_molecular = m2efm_molecular,
    cox_molecular = cox_molecular,
    m2efm_integrated_fit = m2efm_integrated_fit,
    cox_integrated_fit = cox_integrated_fit,
    m2efm_final_risk = m2efm_final_risk,
    cox_final_risk = cox_final_risk
  )

  saveRDS(
    models,
    model_checkpoint,
    compress = "xz"
  )

  cat(
    "Saved full-model checkpoint.\n"
  )
}

m2efm_fit <- models$m2efm_fit
cox_fit <- models$cox_fit
m2efm_molecular <- models$m2efm_molecular
cox_molecular <- models$cox_molecular
m2efm_integrated_fit <- models$m2efm_integrated_fit
cox_integrated_fit <- models$cox_integrated_fit
m2efm_final_risk <- models$m2efm_final_risk
cox_final_risk <- models$cox_final_risk

risk_groups <- setNames(
  lapply(
    cohorts,
    function(cohort) {
      assign_risk_group(
        m2efm_final_risk[[cohort]]
      )
    }
  ),
  cohorts
)


# =====================================================================
# STEP 14 — CALIBRATION CURVES AND IBS
# =====================================================================

cat(
  "\nSTEP 14: calibration and Brier/IBS...\n"
)

step14_result <- try(
  {

    calibration_rows <- list()
    calibration_index <- 1L

    brier_rows <- list()
    brier_index <- 1L

    ibs_rows <- list()
    ibs_index <- 1L

    model_definitions <- list(
      "M2EFM Exp+Clin" = list(
        fit = m2efm_integrated_fit,
        molecular = m2efm_molecular
      ),
      "Cox Exp+Clin" = list(
        fit = cox_integrated_fit,
        molecular = cox_molecular
      )
    )

    five_year_predictions <- list()

    for (model_name in names(
      model_definitions
    )) {

      model_definition <- model_definitions[[model_name]]

      five_year_predictions[[model_name]] <- list()

      for (cohort in cohorts) {

        predicted_at_five <- predict_event_probability(
          model_definition$fit,
          model_definition$molecular[[cohort]],
          clinical[[cohort]],
          times = five_year_days
        )[
          ,
          1L
        ]

        five_year_predictions[[model_name]][[cohort]] <-
          predicted_at_five

        calibration_group <- safe_quantile_groups(
          predicted_at_five,
          groups = 4L
        )

        for (group_level in levels(
          calibration_group
        )) {

          index <- calibration_group ==
            group_level

          observed <- observed_event_probability(
            clinical[[cohort]][[time_col]][index],
            clinical[[cohort]][[event_col]][index],
            five_year_days
          )

          calibration_rows[[calibration_index]] <- data.frame(
            model = model_name,
            cohort = cohort,
            calibration_group = as.integer(
              group_level
            ),
            samples = sum(
              index
            ),
            events = sum(
              clinical[[cohort]][[event_col]][index] ==
                1
            ),
            mean_predicted_event_probability = mean(
              predicted_at_five[index]
            ),
            observed_event_probability = observed,
            calibration_time_days = five_year_days,
            stringsAsFactors = FALSE
          )

          calibration_index <- calibration_index + 1L
        }

        maximum_time <- min(
          ten_year_days,
          max(
            clinical[[cohort]][[time_col]]
          )
        )

        evaluation_times <- unique(
          sort(
            c(
              seq(
                30,
                maximum_time,
                length.out = 80L
              ),
              five_year_days[
                five_year_days <=
                  maximum_time
              ]
            )
          )
        )

        prediction_matrix <- predict_event_probability(
          model_definition$fit,
          model_definition$molecular[[cohort]],
          clinical[[cohort]],
          times = evaluation_times
        )

        brier_values <- vapply(
          seq_along(
            evaluation_times
          ),
          function(j) {
            ipcw_brier(
              time = clinical[[cohort]][[time_col]],
              event = clinical[[cohort]][[event_col]],
              predicted_event_probability = prediction_matrix[
                ,
                j
              ],
              evaluation_time = evaluation_times[[j]]
            )
          },
          numeric(1)
        )

        for (j in seq_along(
          evaluation_times
        )) {

          brier_rows[[brier_index]] <- data.frame(
            model = model_name,
            cohort = cohort,
            time_days = evaluation_times[[j]],
            brier_score = brier_values[[j]],
            stringsAsFactors = FALSE
          )

          brier_index <- brier_index + 1L
        }

        ibs_value <- trapezoid_integral(
          evaluation_times,
          brier_values
        ) /
          (
            max(
              evaluation_times
            ) -
              min(
                evaluation_times
              )
          )

        ibs_rows[[ibs_index]] <- data.frame(
          model = model_name,
          cohort = cohort,
          integration_start_days = min(
            evaluation_times
          ),
          integration_end_days = max(
            evaluation_times
          ),
          integrated_brier_score = ibs_value,
          method = "IPCW censoring-KM; trapezoidal integration",
          stringsAsFactors = FALSE
        )

        ibs_index <- ibs_index + 1L
      }
    }

    calibration_table <- do.call(
      rbind,
      calibration_rows
    )

    brier_table <- do.call(
      rbind,
      brier_rows
    )

    ibs_table <- do.call(
      rbind,
      ibs_rows
    )

    write.csv(
      calibration_table,
      file.path(
        output_dir,
        "14_figureS1_calibration_5year_table.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      brier_table,
      file.path(
        output_dir,
        "14_prediction_error_curves.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      ibs_table,
      file.path(
        output_dir,
        "14_integrated_brier_scores.csv"
      ),
      row.names = FALSE
    )

    plot_calibration <- function(
        file,
        device = c(
          "png",
          "pdf"
        )
    ) {

      device <- match.arg(
        device
      )

      if (device ==
        "png") {

        grDevices::png(
          file,
          width = 2600,
          height = 1250,
          res = 250,
          bg = "white"
        )

      } else {

        grDevices::pdf(
          file,
          width = 10.4,
          height = 5.0
        )
      }

      on.exit(
        grDevices::dev.off(),
        add = TRUE
      )

      old_par <- graphics::par(
        no.readonly = TRUE
      )

      on.exit(
        graphics::par(
          old_par
        ),
        add = TRUE
      )

      graphics::par(
        mfrow = c(
          1,
          2
        ),
        mar = c(
          4.8,
          4.8,
          3.0,
          1.2
        ),
        mgp = c(
          2.8,
          0.8,
          0
        ),
        bty = "l"
      )

      cohort_colors <- c(
        TCGA = "#1B3A8A",
        Kao = "#D95F02",
        Terunuma = "#1B9E77"
      )

      panel_labels <- c(
        "M2EFM Exp+Clin" = "A) M2EFM Exp+Clin",
        "Cox Exp+Clin" = "B) Cox Exp+Clin"
      )

      for (model_name in names(
        model_definitions
      )) {

        current <- calibration_table[
          calibration_table$model ==
            model_name,
          ,
          drop = FALSE
        ]

        plot(
          NA,
          xlim = c(
            0,
            0.75
          ),
          ylim = c(
            0,
            0.75
          ),
          xlab = "Predicted event probability",
          ylab = "Observed event frequency",
          axes = FALSE
        )

        graphics::axis(
          1,
          at = seq(
            0,
            0.75,
            length.out = 5
          ),
          labels = paste0(
            100 *
              seq(
                0,
                0.75,
                length.out = 5
              ),
            " %"
          )
        )

        graphics::axis(
          2,
          at = seq(
            0,
            0.75,
            length.out = 5
          ),
          labels = paste0(
            100 *
              seq(
                0,
                0.75,
                length.out = 5
              ),
            " %"
          ),
          las = 1
        )

        graphics::box(
          bty = "l"
        )

        graphics::abline(
          0,
          1,
          lty = 2,
          lwd = 1.2,
          col = "grey40"
        )

        for (cohort in c(
          "TCGA",
          "Kao",
          "Terunuma"
        )) {

          cohort_data <- current[
            current$cohort ==
              cohort,
            ,
            drop = FALSE
          ]

          cohort_data <- cohort_data[
            order(
              cohort_data$mean_predicted_event_probability
            ),
            ,
            drop = FALSE
          ]

          graphics::lines(
            cohort_data$mean_predicted_event_probability,
            cohort_data$observed_event_probability,
            type = "b",
            lwd = 1.8,
            pch = 16,
            col = cohort_colors[[cohort]]
          )
        }

        graphics::title(
          main = panel_labels[[model_name]],
          adj = 0,
          font.main = 1
        )

        if (model_name ==
          "Cox Exp+Clin") {

          graphics::legend(
            "topleft",
            legend = c(
              "TCGA",
              "Kao",
              "Terunuma"
            ),
            col = cohort_colors[
              c(
                "TCGA",
                "Kao",
                "Terunuma"
              )
            ],
            lty = 1,
            lwd = 1.8,
            pch = 16,
            bty = "n"
          )
        }
      }
    }

    plot_calibration(
      file.path(
        output_dir,
        "14_figureS1_calibration_5year.png"
      ),
      "png"
    )

    plot_calibration(
      file.path(
        output_dir,
        "14_figureS1_calibration_5year.pdf"
      ),
      "pdf"
    )

    list(
      calibration = calibration_table,
      brier = brier_table,
      ibs = ibs_table
    )
  },
  silent = TRUE
)

if (inherits(
  step14_result,
  "try-error"
)) {

  record_status(
    "14",
    "Calibration and IBS",
    "FAILED",
    as.character(
      step14_result
    )
  )

  warning(
    "Step 14 failed: ",
    as.character(
      step14_result
    )
  )

} else {

  record_status(
    "14",
    "Calibration and IBS",
    "COMPLETED",
    "5-year calibration and IPCW IBS written."
  )
}


# =====================================================================
# STEP 15 — KAO SUBGROUP TABLE 2
# =====================================================================

cat(
  "\nSTEP 15: Kao receptor and subgroup C-indices...\n"
)

step15_result <- try(
  {

    kao_clinical_raw <- ready$Kao$clinical
    kao_column_names <- colnames(
      kao_clinical_raw
    )

    er_column <- find_first_column(
      kao_column_names,
      exact_candidates = c(
        "estrogen.receptor.status",
        "er.status",
        "ER.Status",
        "ER"
      ),
      regex_candidates = c(
        "^er[._ ]?status$",
        "estrogen.*receptor"
      )
    )

    pr_column <- find_first_column(
      kao_column_names,
      exact_candidates = c(
        "progesterone.receptor.status",
        "pr.status",
        "PR.Status",
        "PR"
      ),
      regex_candidates = c(
        "^pr[._ ]?status$",
        "progesterone.*receptor"
      )
    )

    her2_column <- find_first_column(
      kao_column_names,
      exact_candidates = c(
        "her2.status",
        "HER2.Status",
        "HER2",
        "her2"
      ),
      regex_candidates = c(
        "^her2[._ ]?status$",
        "her2.*receptor"
      )
    )

    receptor_source <- "clinical_annotations"

    use_clinical <- (
      !is.na(
        er_column
      ) &&
        !is.na(
          pr_column
        ) &&
        !is.na(
          her2_column
        )
    )

    if (use_clinical) {

      er_status <- normalize_binary_status(
        kao_clinical_raw[[er_column]]
      )

      pr_status <- normalize_binary_status(
        kao_clinical_raw[[pr_column]]
      )

      her2_status <- normalize_binary_status(
        kao_clinical_raw[[her2_column]]
      )

      if (
        mean(
          complete.cases(
            er_status,
            pr_status,
            her2_status
          )
        ) <
          0.80
      ) {
        use_clinical <- FALSE
      }
    }

    inference_audit <- NULL

    if (!use_clinical) {

      receptor_source <-
        "ESR1_PGR_ERBB2_two_component_inference"

      required_receptor_genes <- c(
        "ESR1",
        "PGR",
        "ERBB2"
      )

      missing_receptor_genes <- setdiff(
        required_receptor_genes,
        rownames(
          expression$Kao
        )
      )

      if (length(
        missing_receptor_genes
      ) > 0L) {
        stop(
          "Cannot infer Kao receptor status; missing genes: ",
          paste(
            missing_receptor_genes,
            collapse = ", "
          )
        )
      }

      er_inference <- infer_two_component_status(
        expression$Kao[
          "ESR1",
          ,
          drop = TRUE
        ],
        colnames(
          expression$Kao
        ),
        "ESR1_ER"
      )

      pr_inference <- infer_two_component_status(
        expression$Kao[
          "PGR",
          ,
          drop = TRUE
        ],
        colnames(
          expression$Kao
        ),
        "PGR_PR"
      )

      her2_inference <- infer_two_component_status(
        expression$Kao[
          "ERBB2",
          ,
          drop = TRUE
        ],
        colnames(
          expression$Kao
        ),
        "ERBB2_HER2"
      )

      er_status <- er_inference$status
      pr_status <- pr_inference$status
      her2_status <- her2_inference$status

      inference_audit <- rbind(
        er_inference,
        pr_inference,
        her2_inference
      )

      write.csv(
        inference_audit,
        file.path(
          output_dir,
          "15_Kao_receptor_inference_audit.csv"
        ),
        row.names = FALSE
      )
    }

    receptor_table <- data.frame(
      sample_id = rownames(
        clinical$Kao
      ),
      ER_status = er_status,
      PR_status = pr_status,
      HER2_status = her2_status,
      triple_negative = ifelse(
        er_status ==
          "Negative" &
          pr_status ==
            "Negative" &
          her2_status ==
            "Negative",
        "Yes",
        "No"
      ),
      receptor_source = receptor_source,
      stringsAsFactors = FALSE
    )

    write.csv(
      receptor_table,
      file.path(
        output_dir,
        "15_Kao_receptor_status_assignments.csv"
      ),
      row.names = FALSE
    )

    subgroup_definitions <- list(
      Overall = rep(
        TRUE,
        nrow(
          receptor_table
        )
      ),
      `ER-Positive` = receptor_table$ER_status ==
        "Positive",
      `ER-Negative` = receptor_table$ER_status ==
        "Negative",
      `Triple Negative` = receptor_table$triple_negative ==
        "Yes"
    )

    subgroup_rows <- lapply(
      names(
        subgroup_definitions
      ),
      function(subgroup_name) {

        index <- subgroup_definitions[[subgroup_name]]

        index[is.na(
          index
        )] <- FALSE

        cindex_table_row(
          score = m2efm_final_risk$Kao[index],
          clinical = clinical$Kao[
            index,
            ,
            drop = FALSE
          ],
          cohort = "Kao",
          subgroup = subgroup_name
        )
      }
    )

    table2 <- do.call(
      rbind,
      subgroup_rows
    )

    table2$receptor_source <- receptor_source

    table2$interpretation <- if (
      receptor_source ==
        "clinical_annotations"
    ) {
      "Direct available annotations"
    } else {
      "Methodological approximation; exact author probe set unavailable"
    }

    write.csv(
      table2,
      file.path(
        output_dir,
        "15_table2_Kao_subgroup_Cindices.csv"
      ),
      row.names = FALSE
    )

    list(
      table2 = table2,
      receptor_table = receptor_table,
      source = receptor_source
    )
  },
  silent = TRUE
)

if (inherits(
  step15_result,
  "try-error"
)) {

  record_status(
    "15",
    "Kao subgroup Table 2",
    "FAILED",
    as.character(
      step15_result
    )
  )

  warning(
    "Step 15 failed: ",
    as.character(
      step15_result
    )
  )

} else {

  record_status(
    "15",
    "Kao subgroup Table 2",
    "COMPLETED",
    paste(
      "Receptor source:",
      step15_result$source
    )
  )
}


# =====================================================================
# STEP 16 — GO-BP ENRICHMENT AND FIGURE 6 HEATMAP
# =====================================================================

cat(
  "\nSTEP 16: GO-BP enrichment and pathway heatmap...\n"
)

step16_result <- try(
  {

    go_packages <- c(
      "AnnotationDbi",
      "org.Hs.eg.db",
      "GO.db"
    )

    go_ok <- vapply(
      go_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )

    if (!all(
      go_ok
    )) {
      stop(
        "Missing GO packages: ",
        paste(
          names(go_ok)[!go_ok],
          collapse = ", "
        )
      )
    }

    go_mapping <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = common_genes,
      columns = c(
        "GO",
        "ONTOLOGY"
      ),
      keytype = "SYMBOL"
    )

    go_mapping <- go_mapping[
      !is.na(
        go_mapping$GO
      ) &
        go_mapping$ONTOLOGY ==
          "BP",
      c(
        "SYMBOL",
        "GO"
      ),
      drop = FALSE
    ]

    go_mapping <- unique(
      go_mapping
    )

    universe_size <- length(
      common_genes
    )

    signature_in_universe <- intersect(
      signature_genes,
      common_genes
    )

    signature_size <- length(
      signature_in_universe
    )

    term_gene_lists <- split(
      go_mapping$SYMBOL,
      go_mapping$GO
    )

    enrichment_rows <- lapply(
      names(
        term_gene_lists
      ),
      function(go_id) {

        term_genes <- unique(
          intersect(
            term_gene_lists[[go_id]],
            common_genes
          )
        )

        term_size <- length(
          term_genes
        )

        overlap_genes <- intersect(
          signature_in_universe,
          term_genes
        )

        overlap_size <- length(
          overlap_genes
        )

        if (
          term_size < 10L ||
            term_size > 500L ||
            overlap_size < 2L
        ) {
          return(
            NULL
          )
        }

        p_value <- stats::phyper(
          overlap_size -
            1L,
          term_size,
          universe_size -
            term_size,
          signature_size,
          lower.tail = FALSE
        )

        data.frame(
          GO_ID = go_id,
          universe_genes = universe_size,
          term_genes = term_size,
          signature_genes = signature_size,
          overlap_genes = overlap_size,
          p_value = p_value,
          overlap_gene_symbols = paste(
            sort(
              overlap_genes
            ),
            collapse = "/"
          ),
          stringsAsFactors = FALSE
        )
      }
    )

    enrichment_rows <- Filter(
      Negate(
        is.null
      ),
      enrichment_rows
    )

    if (length(
      enrichment_rows
    ) == 0L) {
      stop(
        "No GO-BP terms met minimum enrichment criteria."
      )
    }

    enrichment <- do.call(
      rbind,
      enrichment_rows
    )

    enrichment$p_adjust_BH <- p.adjust(
      enrichment$p_value,
      method = "BH"
    )

    term_annotations <- AnnotationDbi::select(
      GO.db::GO.db,
      keys = unique(
        enrichment$GO_ID
      ),
      columns = "TERM",
      keytype = "GOID"
    )

    term_annotations <- term_annotations[
      !duplicated(
        term_annotations$GOID
      ),
      ,
      drop = FALSE
    ]

    enrichment <- merge(
      enrichment,
      term_annotations[
        ,
        c(
          "GOID",
          "TERM"
        ),
        drop = FALSE
      ],
      by.x = "GO_ID",
      by.y = "GOID",
      all.x = TRUE,
      sort = FALSE
    )

    enrichment <- enrichment[
      order(
        enrichment$p_adjust_BH,
        enrichment$p_value,
        -enrichment$overlap_genes
      ),
      ,
      drop = FALSE
    ]

    write.csv(
      enrichment,
      file.path(
        output_dir,
        "16_Figure6_GO_BP_enrichment.csv"
      ),
      row.names = FALSE
    )

    selected <- head(
      enrichment,
      30L
    )

    selected <- selected[
      !is.na(
        selected$TERM
      ),
      ,
      drop = FALSE
    ]

    if (nrow(
      selected
    ) < 5L) {
      stop(
        "Fewer than five GO-BP terms were available for the heatmap."
      )
    }

    combined_expression <- do.call(
      cbind,
      expression
    )

    combined_expression <- combined_expression[
      signature_in_universe,
      ,
      drop = FALSE
    ]

    gene_means <- matrixStats::rowMeans2(
      combined_expression
    )

    gene_sds <- matrixStats::rowSds(
      combined_expression
    )

    gene_sds[
      !is.finite(
        gene_sds
      ) |
        gene_sds ==
          0
    ] <- 1

    relative_expression <- sweep(
      combined_expression,
      1L,
      gene_means,
      FUN = "-"
    )

    relative_expression <- sweep(
      relative_expression,
      1L,
      gene_sds,
      FUN = "/"
    )

    group_order <- c(
      "TCGA Low",
      "Terunuma Low",
      "Kao Low",
      "TCGA Medium",
      "Terunuma Medium",
      "Kao Medium",
      "TCGA High",
      "Terunuma High",
      "Kao High"
    )

    heatmap_matrix <- matrix(
      NA_real_,
      nrow = length(
        group_order
      ),
      ncol = nrow(
        selected
      ),
      dimnames = list(
        group_order,
        selected$TERM
      )
    )

    for (term_index in seq_len(
      nrow(
        selected
      )
    )) {

      pathway_genes <- strsplit(
        selected$overlap_gene_symbols[[term_index]],
        "/",
        fixed = TRUE
      )[[1L]]

      pathway_genes <- intersect(
        pathway_genes,
        rownames(
          relative_expression
        )
      )

      if (length(
        pathway_genes
      ) == 0L) {
        next
      }

      for (cohort in cohorts) {

        cohort_samples <- colnames(
          expression[[cohort]]
        )

        for (risk_group in risk_levels) {

          selected_samples <- cohort_samples[
            risk_groups[[cohort]] ==
              risk_group
          ]

          row_name <- paste(
            cohort,
            risk_group
          )

          heatmap_matrix[
            row_name,
            term_index
          ] <- mean(
            relative_expression[
              pathway_genes,
              selected_samples,
              drop = FALSE
            ]
          )
        }
      }
    }

    write.csv(
      data.frame(
        group = rownames(
          heatmap_matrix
        ),
        heatmap_matrix,
        check.names = FALSE
      ),
      file.path(
        output_dir,
        "16_Figure6_pathway_mean_relative_expression.csv"
      ),
      row.names = FALSE
    )

    draw_heatmap <- function(
        file,
        device = c(
          "png",
          "pdf"
        )
    ) {

      device <- match.arg(
        device
      )

      if (device ==
        "png") {

        grDevices::png(
          file,
          width = 4200,
          height = 1900,
          res = 260,
          bg = "white"
        )

      } else {

        grDevices::pdf(
          file,
          width = 16.2,
          height = 7.4
        )
      }

      on.exit(
        grDevices::dev.off(),
        add = TRUE
      )

      old_par <- graphics::par(
        no.readonly = TRUE
      )

      on.exit(
        graphics::par(
          old_par
        ),
        add = TRUE
      )

      color_function <- grDevices::colorRampPalette(
        c(
          "#08306B",
          "#6BAED6",
          "#F7FBFF",
          "#FC9272",
          "#99000D"
        )
      )

      graphics::heatmap(
        heatmap_matrix,
        Rowv = NA,
        Colv = NA,
        scale = "none",
        col = color_function(
          120
        ),
        margins = c(
          18,
          13
        ),
        labRow = rownames(
          heatmap_matrix
        ),
        labCol = colnames(
          heatmap_matrix
        ),
        cexRow = 0.9,
        cexCol = 0.55,
        keep.dendro = FALSE
      )

      graphics::title(
        main = "Mean relative expression of enriched GO biological processes"
      )
    }

    draw_heatmap(
      file.path(
        output_dir,
        "16_Figure6_pathway_heatmap.png"
      ),
      "png"
    )

    draw_heatmap(
      file.path(
        output_dir,
        "16_Figure6_pathway_heatmap.pdf"
      ),
      "pdf"
    )

    list(
      enrichment = enrichment,
      selected = selected,
      heatmap = heatmap_matrix
    )
  },
  silent = TRUE
)

if (inherits(
  step16_result,
  "try-error"
)) {

  record_status(
    "16",
    "GO-BP pathway heatmap",
    "FAILED",
    as.character(
      step16_result
    )
  )

  warning(
    "Step 16 failed: ",
    as.character(
      step16_result
    )
  )

} else {

  record_status(
    "16",
    "GO-BP pathway heatmap",
    "COMPLETED",
    paste(
      nrow(
        step16_result$selected
      ),
      "GO-BP terms plotted."
    )
  )
}


# =====================================================================
# STEP 17 — AUTHOR-STYLE OS TABLES
# =====================================================================

cat(
  "\nSTEP 17: assembling OS tables...\n"
)

step17_result <- try(
  {

    table_dir <- file.path(
      output_dir,
      "17_OS_tables"
    )

    dir.create(
      table_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    # Table 1: risk-group hazard ratios.
    table1_source <- file.path(
      results_dir,
      "13C_figure3_hazard_ratios.csv"
    )

    if (file.exists(
      table1_source
    )) {

      table1 <- read.csv(
        table1_source,
        check.names = FALSE
      )

    } else {

      risk_rows <- list()
      risk_index <- 1L

      for (cohort in cohorts) {

        risk_data <- data.frame(
          time = clinical[[cohort]][[time_col]],
          event = clinical[[cohort]][[event_col]],
          risk_group = risk_groups[[cohort]]
        )

        risk_fit <- survival::coxph(
          survival::Surv(
            time,
            event
          ) ~ risk_group,
          data = risk_data,
          ties = "breslow"
        )

        formatted <- format_cox_table(
          risk_fit
        )

        formatted$cohort <- cohort

        risk_rows[[risk_index]] <- formatted

        risk_index <- risk_index + 1L
      }

      table1 <- do.call(
        rbind,
        risk_rows
      )
    }

    write.csv(
      table1,
      file.path(
        table_dir,
        "Table1_OS_risk_group_hazard_ratios.csv"
      ),
      row.names = FALSE
    )

    # Table 2.
    if (
      !inherits(
        step15_result,
        "try-error"
      )
    ) {

      write.csv(
        step15_result$table2,
        file.path(
          table_dir,
          "Table2_Kao_ER_TN_subgroup_Cindices.csv"
        ),
        row.names = FALSE
      )
    }

    # Tables S1/S2: cohort characteristics.
    cohort_characteristics <- function(
        clinical_table,
        cohort
    ) {

      rows <- list()
      index <- 1L
      sample_count <- nrow(
        clinical_table
      )

      rows[[index]] <- data.frame(
        cohort = cohort,
        variable = "Total samples",
        level = "All",
        n = sample_count,
        percent = 100,
        summary_value = as.character(
          sample_count
        ),
        stringsAsFactors = FALSE
      )

      index <- index + 1L

      deaths <- sum(
        clinical_table[[event_col]] ==
          1
      )

      rows[[index]] <- data.frame(
        cohort = cohort,
        variable = "Deaths",
        level = "Event",
        n = deaths,
        percent = 100 *
          deaths /
          sample_count,
        summary_value = as.character(
          deaths
        ),
        stringsAsFactors = FALSE
      )

      index <- index + 1L

      rows[[index]] <- data.frame(
        cohort = cohort,
        variable = "Age at diagnosis",
        level = "Mean (SD)",
        n = sample_count,
        percent = NA_real_,
        summary_value = sprintf(
          "%.2f (%.2f)",
          mean(
            clinical_table[[age_col]]
          ),
          sd(
            clinical_table[[age_col]]
          )
        ),
        stringsAsFactors = FALSE
      )

      index <- index + 1L

      stage_counts <- table(
        clinical_table[[stage_col]],
        useNA = "ifany"
      )

      for (stage_name in names(
        stage_counts
      )) {

        count <- as.integer(
          stage_counts[[stage_name]]
        )

        rows[[index]] <- data.frame(
          cohort = cohort,
          variable = "Pathologic stage",
          level = stage_name,
          n = count,
          percent = 100 *
            count /
            sample_count,
          summary_value = as.character(
            count
          ),
          stringsAsFactors = FALSE
        )

        index <- index + 1L
      }

      do.call(
        rbind,
        rows
      )
    }

    characteristics <- do.call(
      rbind,
      lapply(
        cohorts,
        function(cohort) {
          cohort_characteristics(
            clinical[[cohort]],
            cohort
          )
        }
      )
    )

    write.csv(
      characteristics[
        characteristics$cohort ==
          "TCGA",
        ,
        drop = FALSE
      ],
      file.path(
        table_dir,
        "TableS1_TCGA_characteristics_replication.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      characteristics[
        characteristics$cohort %in%
          c(
            "Terunuma",
            "Kao"
          ),
        ,
        drop = FALSE
      ],
      file.path(
        table_dir,
        "TableS2_external_characteristics_replication.csv"
      ),
      row.names = FALSE
    )

    # Table S4: available expression signature.
    table_s4 <- data.frame(
      feature = signature_genes,
      feature_type = "external-common M2EFM expression gene",
      available_in_TCGA = signature_genes %in%
        rownames(
          expression$TCGA
        ),
      available_in_Terunuma = signature_genes %in%
        rownames(
          expression$Terunuma
        ),
      available_in_Kao = signature_genes %in%
        rownames(
          expression$Kao
        ),
      stringsAsFactors = FALSE
    )

    write.csv(
      table_s4,
      file.path(
        table_dir,
        "TableS4_available_M2EFM_expression_features.csv"
      ),
      row.names = FALSE
    )

    # Table S5: prior Figure 2A summary, discovered by filename.
    figure2a_candidates <- list.files(
      results_dir,
      pattern = "^08.*summary.*\\.csv$",
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(
      figure2a_candidates
    ) > 0L) {

      table_s5 <- read.csv(
        sort(
          figure2a_candidates
        )[[1L]],
        check.names = FALSE
      )

      write.csv(
        table_s5,
        file.path(
          table_dir,
          "TableS5_meta_dimensional_OS_models_replication.csv"
        ),
        row.names = FALSE
      )

      record_status(
        "17",
        "Table S5",
        "COMPLETED",
        basename(
          sort(
            figure2a_candidates
          )[[1L]]
        )
      )

    } else {

      record_status(
        "17",
        "Table S5",
        "SKIPPED",
        "No Step 08 summary CSV discovered."
      )
    }

    # Table S6: prior Figure 2B summary.
    table_s6_source <- file.path(
      results_dir,
      "11_figure2B_summary_100splits.csv"
    )

    if (file.exists(
      table_s6_source
    )) {

      table_s6 <- read.csv(
        table_s6_source,
        check.names = FALSE
      )

      write.csv(
        table_s6,
        file.path(
          table_dir,
          "TableS6_expression_only_OS_models_7_available.csv"
        ),
        row.names = FALSE
      )

    } else {

      record_status(
        "17",
        "Table S6",
        "SKIPPED",
        "Step 11 summary file missing."
      )
    }

    # Table S7: final full-TCGA M2EFM integrated model HRs.
    table_s7 <- format_cox_table(
      m2efm_integrated_fit
    )

    write.csv(
      table_s7,
      file.path(
        table_dir,
        "TableS7_final_M2EFM_OS_model_hazard_ratios.csv"
      ),
      row.names = FALSE
    )

    # Final-model C-index table.
    final_cindex_rows <- do.call(
      rbind,
      lapply(
        cohorts,
        function(cohort) {
          cindex_table_row(
            m2efm_final_risk[[cohort]],
            clinical[[cohort]],
            cohort,
            "Overall"
          )
        }
      )
    )

    write.csv(
      final_cindex_rows,
      file.path(
        table_dir,
        "Final_full_TCGA_model_Cindices.csv"
      ),
      row.names = FALSE
    )

    # Optional XLSX workbook.
    workbook_path <- file.path(
      output_dir,
      "17_all_OS_tables.xlsx"
    )

    if (requireNamespace(
      "openxlsx",
      quietly = TRUE
    )) {

      workbook <- openxlsx::createWorkbook()

      sheets <- list(
        Table1 = table1,
        TableS1 = characteristics[
          characteristics$cohort ==
            "TCGA",
          ,
          drop = FALSE
        ],
        TableS2 = characteristics[
          characteristics$cohort %in%
            c(
              "Terunuma",
              "Kao"
            ),
          ,
          drop = FALSE
        ],
        TableS4 = table_s4,
        TableS7 = table_s7,
        Final_Cindex = final_cindex_rows
      )

      if (
        !inherits(
          step15_result,
          "try-error"
        )
      ) {
        sheets$Table2 <- step15_result$table2
      }

      if (exists(
        "table_s5"
      )) {
        sheets$TableS5 <- table_s5
      }

      if (exists(
        "table_s6"
      )) {
        sheets$TableS6 <- table_s6
      }

      if (
        !inherits(
          step14_result,
          "try-error"
        )
      ) {
        sheets$Calibration <- step14_result$calibration
        sheets$IBS <- step14_result$ibs
      }

      for (sheet_name in names(
        sheets
      )) {

        openxlsx::addWorksheet(
          workbook,
          sheet_name
        )

        openxlsx::writeData(
          workbook,
          sheet_name,
          sheets[[sheet_name]]
        )

        openxlsx::freezePane(
          workbook,
          sheet_name,
          firstRow = TRUE
        )

        openxlsx::setColWidths(
          workbook,
          sheet_name,
          cols = seq_len(
            ncol(
              sheets[[sheet_name]]
            )
          ),
          widths = "auto"
        )
      }

      openxlsx::saveWorkbook(
        workbook,
        workbook_path,
        overwrite = TRUE
      )

      workbook_note <- basename(
        workbook_path
      )

    } else {

      workbook_note <-
        "openxlsx not installed; CSV tables were still written."
    }

    list(
      table_dir = table_dir,
      workbook_note = workbook_note,
      final_cindices = final_cindex_rows,
      table_s7 = table_s7
    )
  },
  silent = TRUE
)

if (inherits(
  step17_result,
  "try-error"
)) {

  record_status(
    "17",
    "OS table assembly",
    "FAILED",
    as.character(
      step17_result
    )
  )

  warning(
    "Step 17 failed: ",
    as.character(
      step17_result
    )
  )

} else {

  record_status(
    "17",
    "OS table assembly",
    "COMPLETED",
    step17_result$workbook_note
  )
}


# =====================================================================
# FINAL STATUS, AUDIT, AND SESSION INFORMATION
# =====================================================================

status_table <- do.call(
  rbind,
  status_rows
)

write.csv(
  status_table,
  file.path(
    output_dir,
    "14_17_completion_status.csv"
  ),
  row.names = FALSE
)

model_audit <- data.frame(
  metric = c(
    "input_file",
    "input_md5",
    "TCGA_samples",
    "Terunuma_samples",
    "Kao_samples",
    "common_expression_genes",
    "M2EFM_expression_genes",
    "M2EFM_lambda_min",
    "Cox_lambda_min",
    "calibration_time_days",
    "IBS_upper_limit_days",
    "Kao_receptor_exact_probe_set_available",
    "Figure6_database"
  ),
  value = c(
    normalizePath(
      input_file,
      winslash = "/",
      mustWork = TRUE
    ),
    input_md5,
    nrow(
      clinical$TCGA
    ),
    nrow(
      clinical$Terunuma
    ),
    nrow(
      clinical$Kao
    ),
    length(
      common_genes
    ),
    length(
      signature_genes
    ),
    m2efm_fit$lambda.min,
    cox_fit$lambda.min,
    five_year_days,
    ten_year_days,
    FALSE,
    paste0(
      "Installed GO.db ",
      if (
        requireNamespace(
          "GO.db",
          quietly = TRUE
        )
      ) {
        as.character(
          utils::packageVersion(
            "GO.db"
          )
        )
      } else {
        "not available"
      }
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  model_audit,
  file.path(
    output_dir,
    "14_17_model_audit.csv"
  ),
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

capture.output(
  {
    cat(
      "STEPS 14–17 — REMAINING OS ANALYSES\n\n"
    )

    cat(
      "COMPLETION STATUS\n"
    )

    print(
      status_table
    )

    cat(
      "\nMODEL AUDIT\n"
    )

    print(
      model_audit
    )

    if (
      !inherits(
        step14_result,
        "try-error"
      )
    ) {

      cat(
        "\nINTEGRATED BRIER SCORES\n"
      )

      print(
        step14_result$ibs
      )
    }

    if (
      !inherits(
        step15_result,
        "try-error"
      )
    ) {

      cat(
        "\nKAO TABLE 2\n"
      )

      print(
        step15_result$table2
      )
    }

    if (
      !inherits(
        step17_result,
        "try-error"
      )
    ) {

      cat(
        "\nFINAL MODEL C-INDICES\n"
      )

      print(
        step17_result$final_cindices
      )

      cat(
        "\nFINAL MODEL HAZARD RATIOS\n"
      )

      print(
        step17_result$table_s7
      )
    }

    cat(
      "\nElapsed minutes: ",
      elapsed_minutes,
      "\n",
      sep = ""
    )

    cat(
      "\nSESSION INFO\n"
    )

    sessionInfo()
  },
  file = file.path(
    output_dir,
    "14_17_sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n"
)

cat(
  "STEPS 14–17 FINISHED\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nCOMPLETION STATUS\n"
)

print(
  status_table,
  row.names = FALSE
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n",
  sep = ""
)

cat(
  "\nElapsed minutes: ",
  round(
    elapsed_minutes,
    3
  ),
  "\n",
  sep = ""
)

cat(
  "\nFinished.\n"
)
