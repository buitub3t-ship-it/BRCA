# ============================================================
# 03_train_internal_9models.R
#
# Train the 9 models shown in Figure 2A:
#   1. Clin
#   2. M2EFM Meth+Exp+Clin
#   3. M2EFM Exp+Clin
#   4. Cox Meth+Exp+Clin
#   5. rorS+Clin
#   6. M2EFM Meth+Exp
#   7. M2EFM Exp
#   8. Cox Meth+Exp
#   9. rorS
#
# This script uses the cache objects already produced by scripts 01 and 02.
# Run from the M2EFM_OS_TCGA_KAO_TERU project root.
# ============================================================

source("config.R")
source(file.path("R", "io_and_preprocessing.R"))
source(file.path("R", "m2eqtl.R"))
source(file.path("R", "modeling.R"))

required_packages <- c(
  "data.table", "survival", "survcomp",
  "glmnet", "matrixStats"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

needed_cache <- c(
  "tcga_exp_common.rds",
  "tcga_meth_beta.rds",
  "tcga_clin_clean.rds",
  "m2eqtl_discovery.rds"
)

missing_cache <- needed_cache[
  !file.exists(file.path(CACHE_DIR, needed_cache))
]

if (length(missing_cache) > 0L) {
  stop(
    "Missing cache files:\n- ",
    paste(missing_cache, collapse = "\n- "),
    "\nRun scripts 01 and 02 first."
  )
}

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

TIME_COL <- "OVERALL.SURVIVAL"
EVENT_COL <- "overall.survival.indicator"
CLINICAL_COVARIATES <- c("pathologic_stage", "age.Dx")

# TRUE: reuse the already completed 3-model result produced by
# 03_train_internal_full.R, rather than fitting those models again.
REUSE_EXISTING_FULL_M2EFM <- TRUE

# TRUE: reuse expression-only M2EFM cache if it already exists.
REUSE_EXISTING_EXP_M2EFM <- TRUE

# TRUE: resume Cox/rorS splits from a partial checkpoint.
RESUME_PARTIAL_COMPARISONS <- TRUE

# The paper selected approximately as many variable methylation probes
# as expression genes for the Cox-Ridge comparison. The current matrix
# contains fewer methylation probes than genes, so all available probes
# may be retained.
COX_METHYLATION_PROBE_LIMIT <- nrow(
  readRDS(file.path(CACHE_DIR, "tcga_exp_common.rds"))
)

OUTPUT_METRICS <- file.path(
  RESULTS_DIR,
  "03_internal_9models_monte_carlo_metrics.csv"
)

OUTPUT_SUMMARY <- file.path(
  RESULTS_DIR,
  "03_internal_9models_summary.csv"
)

OUTPUT_AUDIT <- file.path(
  RESULTS_DIR,
  "03_internal_9models_audit.csv"
)

OUTPUT_RDS <- file.path(
  CACHE_DIR,
  "internal_9models_monte_carlo.rds"
)

PARTIAL_RDS <- file.path(
  CACHE_DIR,
  "internal_9models_comparisons_partial.rds"
)

PARTIAL_ERRORS_RDS <- file.path(
  CACHE_DIR,
  "internal_9models_comparison_errors_partial.rds"
)

OUTPUT_ERRORS <- file.path(
  RESULTS_DIR,
  "03_internal_9models_comparison_errors.csv"
)

EXP_M2EFM_RDS <- file.path(
  CACHE_DIR,
  "internal_exp_monte_carlo.rds"
)

RORS_RDS <- file.path(
  CACHE_DIR,
  "tcga_pam50_rors.rds"
)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

as_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

normalize_name <- function(x) {
  gsub("[^a-z0-9]", "", tolower(as.character(x)))
}

assert_columns <- function(x, columns, object_name) {
  absent <- setdiff(columns, colnames(x))
  if (length(absent) > 0L) {
    stop(
      object_name, " is missing columns: ",
      paste(absent, collapse = ", "),
      "\nAvailable columns:\n",
      paste(colnames(x), collapse = ", ")
    )
  }
}

safe_cindex <- function(score, time, event) {
  score <- as_numeric_safe(score)
  time <- as_numeric_safe(time)
  event <- as_numeric_safe(event)

  keep <- is.finite(score) & is.finite(time) & is.finite(event)

  if (sum(keep) < 3L || length(unique(event[keep])) < 2L) {
    return(NA_real_)
  }

  value <- tryCatch(
    {
      z <- survcomp::concordance.index(
        x = score[keep],
        surv.time = time[keep],
        surv.event = event[keep]
      )
      as.numeric(z$c.index)
    },
    error = function(e) NA_real_
  )

  if (!is.finite(value)) {
    value <- tryCatch(
      {
        fit <- survival::concordance(
          survival::Surv(time[keep], event[keep]) ~ score[keep],
          reverse = TRUE
        )
        as.numeric(fit$concordance)
      },
      error = function(e) NA_real_
    )
  }

  value
}

make_survival <- function(clin) {
  survival::Surv(
    time = as_numeric_safe(clin[[TIME_COL]]),
    event = as_numeric_safe(clin[[EVENT_COL]])
  )
}

clinical_formula <- function(response = "outcome") {
  stats::as.formula(
    paste(
      response,
      "~",
      paste(CLINICAL_COVARIATES, collapse = " + ")
    )
  )
}

combined_formula <- function(response = "outcome") {
  stats::as.formula(
    paste(
      response,
      "~ molecular_risk +",
      paste(CLINICAL_COVARIATES, collapse = " + ")
    )
  )
}

prepare_clinical <- function(clin, sample_ids) {
  clin <- clin[sample_ids, , drop = FALSE]

  assert_columns(
    clin,
    c(TIME_COL, EVENT_COL, CLINICAL_COVARIATES),
    "TCGA clinical"
  )

  clin[[TIME_COL]] <- as_numeric_safe(clin[[TIME_COL]])
  clin[[EVENT_COL]] <- as_numeric_safe(clin[[EVENT_COL]])
  clin[["age.Dx"]] <- as_numeric_safe(clin[["age.Dx"]])

  if (exists("STAGE_LEVELS", inherits = TRUE)) {
    clin[["pathologic_stage"]] <- factor(
      as.character(clin[["pathologic_stage"]]),
      levels = STAGE_LEVELS
    )
  } else {
    clin[["pathologic_stage"]] <- factor(
      as.character(clin[["pathologic_stage"]])
    )
  }

  complete <- stats::complete.cases(
    clin[, c(TIME_COL, EVENT_COL, CLINICAL_COVARIATES), drop = FALSE]
  )

  clin <- clin[complete, , drop = FALSE]

  observed_events <- sort(unique(clin[[EVENT_COL]]))
  if (!all(observed_events %in% c(0, 1))) {
    stop(
      EVENT_COL, " must be coded 0/1. Observed: ",
      paste(observed_events, collapse = ", ")
    )
  }

  clin
}

remove_invalid_features <- function(x, label) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  valid <- apply(
    x,
    1L,
    function(v) {
      all(is.finite(v)) &&
        stats::sd(v) > 0
    }
  )

  if (any(!valid)) {
    message(
      label, ": removing ", sum(!valid),
      " non-finite or zero-variance features."
    )
  }

  x[valid, , drop = FALSE]
}

replay_split_rng <- function(split_id, n_samples, train_size) {
  # Recreate the RNG state used by the author's evaluate():
  # set.seed(i), followed by sample() for the 70/30 split.
  set.seed(split_id)
  invisible(
    sample.int(
      n = n_samples,
      size = train_size,
      replace = FALSE
    )
  )
}

fit_ridge_risk <- function(
    sample_by_feature,
    outcome,
    split_id,
    n_total,
    train_size
) {
  replay_split_rng(split_id, n_total, train_size)

  fit <- glmnet::cv.glmnet(
    x = sample_by_feature,
    y = outcome,
    family = "cox",
    alpha = 0,
    nfolds = GLMNET_CV_FOLDS,
    standardize = FALSE
  )

  list(
    fit = fit,
    risk = as.numeric(
      stats::predict(
        fit,
        newx = sample_by_feature,
        s = "lambda.min",
        type = "response"
      )
    )
  )
}

fit_clinical_model <- function(clin_train) {
  outcome <- make_survival(clin_train)

  survival::coxph(
    formula = clinical_formula("outcome"),
    data = clin_train,
    ties = "efron",
    model = TRUE,
    x = TRUE
  )
}

fit_combined_model <- function(molecular_risk, clin_train) {
  combined <- data.frame(
    molecular_risk = as.numeric(molecular_risk),
    clin_train[, CLINICAL_COVARIATES, drop = FALSE],
    check.names = FALSE
  )

  outcome <- make_survival(clin_train)

  survival::coxph(
    formula = combined_formula("outcome"),
    data = combined,
    ties = "efron",
    model = TRUE,
    x = TRUE
  )
}

predict_clinical_model <- function(fit, clin_test) {
  as.numeric(
    stats::predict(
      fit,
      newdata = clin_test,
      type = "lp"
    )
  )
}

predict_combined_model <- function(
    fit,
    molecular_risk,
    clin_test
) {
  new_data <- data.frame(
    molecular_risk = as.numeric(molecular_risk),
    clin_test[, CLINICAL_COVARIATES, drop = FALSE],
    check.names = FALSE
  )

  as.numeric(
    stats::predict(
      fit,
      newdata = new_data,
      type = "lp"
    )
  )
}

canonicalize_existing_full_metrics <- function(metrics) {
  metrics <- data.table::as.data.table(metrics)

  required <- c("split", "cohort", "model", "c_index")
  absent <- setdiff(required, colnames(metrics))

  if (length(absent) > 0L) {
    stop(
      "Existing full M2EFM metrics are missing columns: ",
      paste(absent, collapse = ", ")
    )
  }

  mapping <- c(
    "Clinical" = "Clin",
    "M2EFM_Molecular" = "M2EFM Meth+Exp",
    "M2EFM_ExpClin" = "M2EFM Meth+Exp+Clin"
  )

  metrics <- metrics[
    cohort == "TCGA_test" &
      model %in% names(mapping)
  ]

  metrics[, model := unname(mapping[model])]
  metrics[, cohort := "TCGA_test"]

  metrics[, .(
    split = as.integer(split),
    cohort,
    model,
    c_index = as.numeric(c_index)
  )]
}

canonicalize_exp_metrics <- function(metrics) {
  metrics <- data.table::as.data.table(metrics)

  mapping <- c(
    "M2EFM_Molecular" = "M2EFM Exp",
    "M2EFM_ExpClin" = "M2EFM Exp+Clin"
  )

  metrics <- metrics[
    cohort == "TCGA_test" &
      model %in% names(mapping)
  ]

  metrics[, model := unname(mapping[model])]
  metrics[, cohort := "TCGA_test"]

  metrics[, .(
    split = as.integer(split),
    cohort,
    model,
    c_index = as.numeric(c_index)
  )]
}

find_numeric_rors_column <- function(clin) {
  keys <- normalize_name(colnames(clin))
  preferred <- normalize_name(c(
    "rorS", "ror_score", "pam50_rors",
    "pam50risk", "pam50_risk_score"
  ))

  hits <- which(keys %in% preferred)

  for (index in hits) {
    values <- as_numeric_safe(clin[[index]])
    if (sum(is.finite(values)) >= 50L) {
      names(values) <- rownames(clin)
      return(values)
    }
  }

  NULL
}

compute_or_load_rors <- function(tcga_exp, tcga_clin) {
  # Reuse a successfully computed score if available.
  if (file.exists(RORS_RDS)) {
    score <- readRDS(RORS_RDS)
    original_names <- names(score)
    score <- as_numeric_safe(score)
    names(score) <- original_names

    if (
      length(score) > 0L &&
        !is.null(names(score)) &&
        any(is.finite(score))
    ) {
      message("Loaded cached PAM50 rorS: ", RORS_RDS)
      return(score)
    }

    warning(
      "Ignoring invalid rorS cache and recomputing: ",
      RORS_RDS
    )
  }

  # Prefer an rorS column already present in the clinical file.
  score <- find_numeric_rors_column(tcga_clin)

  if (!is.null(score)) {
    saveRDS(score, RORS_RDS)
    message("Using an existing numeric rorS column from TCGA clinical data.")
    return(score)
  }

  required_rors_packages <- c(
    "genefu",
    "org.Hs.eg.db",
    "AnnotationDbi"
  )

  missing_rors_packages <- required_rors_packages[
    !vapply(
      required_rors_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_rors_packages) > 0L) {
    stop(
      "Cannot compute PAM50 rorS. Missing packages: ",
      paste(missing_rors_packages, collapse = ", "),
      "\nInstall them, restart R, and rerun:",
      "\nBiocManager::install(c(",
      paste(
        sprintf('"%s"', required_rors_packages),
        collapse = ", "
      ),
      "))"
    )
  }

  message(
    "Computing PAM50 rorS directly with genefu::rorS() ..."
  )

  # This reproduces the rorS part of M2EFM::get_pam50(), but deliberately
  # omits molecular.subtyping(). The 9-model analysis needs only the
  # continuous rorS score, not the PAM50 subtype calls.
  alias_map <- org.Hs.eg.db::org.Hs.egALIAS2EG
  mapped_aliases <- AnnotationDbi::mappedkeys(alias_map)
  alias_to_entrez <- AnnotationDbi::as.list(
    alias_map[mapped_aliases]
  )

  matched_genes <- intersect(
    rownames(tcga_exp),
    names(alias_to_entrez)
  )

  if (length(matched_genes) < 50L) {
    stop(
      "Too few expression genes could be mapped for rorS: ",
      length(matched_genes)
    )
  }

  exp_for_rors <- tcga_exp[
    matched_genes,
    ,
    drop = FALSE
  ]

  entrez_ids <- vapply(
    alias_to_entrez[matched_genes],
    function(z) {
      if (length(z) == 0L || is.na(z[[1L]])) {
        NA_character_
      } else {
        as.character(z[[1L]])
      }
    },
    character(1)
  )

  valid_annotation <- !is.na(entrez_ids) &
    nzchar(entrez_ids)

  exp_for_rors <- exp_for_rors[
    valid_annotation,
    ,
    drop = FALSE
  ]
  entrez_ids <- entrez_ids[valid_annotation]

  annot <- data.frame(
    GeneSymbol = rownames(exp_for_rors),
    EntrezGene.ID = entrez_ids,
    row.names = rownames(exp_for_rors),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rors_result <- genefu::rorS(
    data = t(exp_for_rors),
    annot = annot,
    do.mapping = FALSE,
    verbose = FALSE
  )

  if (
    is.null(rors_result$score) ||
      length(rors_result$score) == 0L
  ) {
    stop("genefu::rorS() returned no continuous score.")
  }

  score_raw <- rors_result$score
  score <- as_numeric_safe(score_raw)
  names(score) <- gsub(
    "\\.",
    "-",
    names(score_raw)
  )

  if (sum(is.finite(score)) < 50L) {
    stop(
      "genefu::rorS() produced too few finite scores: ",
      sum(is.finite(score))
    )
  }

  saveRDS(score, RORS_RDS)

  message(
    "Saved PAM50 rorS for ",
    sum(is.finite(score)),
    " samples: ",
    RORS_RDS
  )

  score
}

align_named_score <- function(score, sample_ids) {
  if (is.null(names(score))) {
    stop("The rorS vector has no sample names.")
  }

  direct <- score[sample_ids]

  if (sum(is.finite(direct)) == length(sample_ids)) {
    return(as.numeric(direct))
  }

  normalized_names <- gsub("\\.", "-", names(score))
  names(score) <- normalized_names
  direct <- score[sample_ids]

  if (sum(is.finite(direct)) < length(sample_ids)) {
    missing_ids <- sample_ids[!is.finite(direct)]
    stop(
      "rorS is missing for ", length(missing_ids),
      " matched TCGA samples. First missing IDs:\n",
      paste(utils::head(missing_ids, 10L), collapse = "\n")
    )
  }

  as.numeric(direct)
}

validate_model_counts <- function(metrics, expected_splits) {
  models <- c(
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

  audit <- data.table::data.table(model = models)
  counts <- metrics[
    is.finite(c_index),
    .(n_valid = .N),
    by = model
  ]

  audit <- merge(
    audit,
    counts,
    by = "model",
    all.x = TRUE,
    sort = FALSE
  )

  audit[is.na(n_valid), n_valid := 0L]
  audit[, expected := expected_splits]
  audit[, complete := n_valid == expected]

  audit
}

# ------------------------------------------------------------
# Load prepared data
# ------------------------------------------------------------

tcga_exp <- readRDS(
  file.path(CACHE_DIR, "tcga_exp_common.rds")
)
tcga_meth_beta <- readRDS(
  file.path(CACHE_DIR, "tcga_meth_beta.rds")
)
tcga_clin <- readRDS(
  file.path(CACHE_DIR, "tcga_clin_clean.rds")
)
m2e <- readRDS(
  file.path(CACHE_DIR, "m2eqtl_discovery.rds")
)

assert_columns(
  tcga_clin,
  c(TIME_COL, EVENT_COL, CLINICAL_COVARIATES),
  "tcga_clin_clean.rds"
)

missing_trans <- setdiff(
  m2e$trans_genes,
  rownames(tcga_exp)
)
missing_probes <- setdiff(
  m2e$methylation_probes,
  rownames(tcga_meth_beta)
)

if (length(missing_trans) > 0L) {
  stop(
    "Selected trans genes missing from tcga_exp_common.rds:\n",
    paste(missing_trans, collapse = ", ")
  )
}

if (length(missing_probes) > 0L) {
  stop(
    "Selected methylation probes missing from tcga_meth_beta.rds:\n",
    paste(missing_probes, collapse = ", ")
  )
}

selected_exp_raw <- tcga_exp[
  m2e$trans_genes,
  ,
  drop = FALSE
]

selected_meth_m <- beta_to_m(
  tcga_meth_beta[
    m2e$methylation_probes,
    ,
    drop = FALSE
  ]
)

selected_exp_scale <- scale_features_reference(
  selected_exp_raw
)
selected_meth_scale <- scale_features_reference(
  selected_meth_m
)

rownames(selected_exp_scale$values) <- paste0(
  "EXP__",
  rownames(selected_exp_scale$values)
)
rownames(selected_meth_scale$values) <- paste0(
  "METH__",
  rownames(selected_meth_scale$values)
)

common_samples <- Reduce(
  intersect,
  list(
    colnames(selected_exp_scale$values),
    colnames(selected_meth_scale$values),
    rownames(tcga_clin)
  )
)

tcga_clin_matched <- prepare_clinical(
  tcga_clin,
  common_samples
)

common_samples <- rownames(tcga_clin_matched)

selected_exp <- selected_exp_scale$values[
  ,
  common_samples,
  drop = FALSE
]

selected_full <- rbind(
  selected_exp_scale$values[
    ,
    common_samples,
    drop = FALSE
  ],
  selected_meth_scale$values[
    ,
    common_samples,
    drop = FALSE
  ]
)

selected_exp <- remove_invalid_features(
  selected_exp,
  "Selected expression"
)
selected_full <- remove_invalid_features(
  selected_full,
  "Selected Meth+Exp"
)

message("Matched TCGA samples: ", length(common_samples))
message("Selected expression features: ", nrow(selected_exp))
message("Selected integrated features: ", nrow(selected_full))

# ------------------------------------------------------------
# 1-3 and 6-7: M2EFM models
# ------------------------------------------------------------

existing_full_csv <- file.path(
  RESULTS_DIR,
  "03_internal_full_monte_carlo_metrics.csv"
)

full_metrics <- NULL

if (
  REUSE_EXISTING_FULL_M2EFM &&
    file.exists(existing_full_csv)
) {
  candidate <- data.table::fread(existing_full_csv)
  candidate <- canonicalize_existing_full_metrics(candidate)

  candidate_audit <- candidate[
    is.finite(c_index),
    .N,
    by = model
  ]

  expected_full <- c(
    "Clin",
    "M2EFM Meth+Exp",
    "M2EFM Meth+Exp+Clin"
  )

  valid_existing <- all(
    expected_full %in% candidate_audit$model
  ) && all(
    candidate_audit[
      model %in% expected_full,
      N
    ] == N_MONTE_CARLO_SPLITS
  )

  if (valid_existing) {
    full_metrics <- candidate
    message(
      "Reusing existing full M2EFM metrics: ",
      existing_full_csv
    )
  }
}

if (is.null(full_metrics)) {
  message("Training full M2EFM models ...")

  full_result <- monte_carlo_m2efm(
    training_feature_by_sample = selected_full,
    training_clinical = tcga_clin_matched,
    validation_sets = list(),
    n_splits = N_MONTE_CARLO_SPLITS
  )

  full_metrics <- canonicalize_existing_full_metrics(
    full_result$metrics
  )

  saveRDS(
    full_result,
    file.path(
      CACHE_DIR,
      "internal_full_monte_carlo_9models_run.rds"
    )
  )
}

exp_result <- NULL

if (
  REUSE_EXISTING_EXP_M2EFM &&
    file.exists(EXP_M2EFM_RDS)
) {
  exp_result <- readRDS(EXP_M2EFM_RDS)

  if (
    is.null(exp_result$metrics) ||
      nrow(exp_result$metrics) == 0L
  ) {
    exp_result <- NULL
  } else {
    message(
      "Reusing expression-only M2EFM cache: ",
      EXP_M2EFM_RDS
    )
  }
}

if (is.null(exp_result)) {
  message("Training expression-only M2EFM models ...")

  exp_result <- monte_carlo_m2efm(
    training_feature_by_sample = selected_exp,
    training_clinical = tcga_clin_matched,
    validation_sets = list(),
    n_splits = N_MONTE_CARLO_SPLITS
  )

  saveRDS(exp_result, EXP_M2EFM_RDS)
}

exp_metrics <- canonicalize_exp_metrics(
  exp_result$metrics
)

# ------------------------------------------------------------
# 4-5 and 8-9: Cox-Ridge and rorS models
# ------------------------------------------------------------

message("Preparing Cox-Ridge all-feature matrix ...")

cox_exp_raw <- tcga_exp[
  ,
  common_samples,
  drop = FALSE
]

cox_exp_scale <- scale_features_reference(
  cox_exp_raw
)
rownames(cox_exp_scale$values) <- paste0(
  "EXPALL__",
  rownames(cox_exp_scale$values)
)

cox_meth_m <- beta_to_m(
  tcga_meth_beta[
    ,
    common_samples,
    drop = FALSE
  ]
)

meth_mad <- matrixStats::rowMads(
  as.matrix(cox_meth_m),
  na.rm = TRUE
)

valid_mad <- is.finite(meth_mad)
meth_order <- order(
  meth_mad[valid_mad],
  decreasing = TRUE
)

valid_probe_names <- rownames(cox_meth_m)[valid_mad]

n_cox_meth <- min(
  length(valid_probe_names),
  as.integer(COX_METHYLATION_PROBE_LIMIT)
)

selected_cox_probes <- valid_probe_names[
  meth_order[seq_len(n_cox_meth)]
]

cox_meth_scale <- scale_features_reference(
  cox_meth_m[
    selected_cox_probes,
    ,
    drop = FALSE
  ]
)
rownames(cox_meth_scale$values) <- paste0(
  "METHALL__",
  rownames(cox_meth_scale$values)
)

cox_full <- rbind(
  cox_exp_scale$values,
  cox_meth_scale$values
)

cox_full <- remove_invalid_features(
  cox_full,
  "Cox Meth+Exp"
)

cox_sample_by_feature <- t(cox_full)
storage.mode(cox_sample_by_feature) <- "double"

message("Cox expression features: ", nrow(cox_exp_scale$values))
message("Cox methylation features: ", nrow(cox_meth_scale$values))
message("Cox total features: ", ncol(cox_sample_by_feature))

rm(
  cox_exp_raw,
  cox_meth_m,
  cox_full
)
invisible(gc())

rors_named <- compute_or_load_rors(
  tcga_exp,
  tcga_clin
)

rors <- align_named_score(
  rors_named,
  common_samples
)
names(rors) <- common_samples

comparison_metrics <- data.table::data.table(
  split = integer(),
  cohort = character(),
  model = character(),
  c_index = numeric()
)

comparison_errors <- data.table::data.table(
  split = integer(),
  model_group = character(),
  error_message = character()
)

if (
  RESUME_PARTIAL_COMPARISONS &&
    file.exists(PARTIAL_RDS)
) {
  partial <- readRDS(PARTIAL_RDS)

  if (is.data.frame(partial) || data.table::is.data.table(partial)) {
    comparison_metrics <- data.table::as.data.table(partial)
    message(
      "Resuming comparison metrics from: ",
      PARTIAL_RDS
    )
  }
}

if (
  RESUME_PARTIAL_COMPARISONS &&
    file.exists(PARTIAL_ERRORS_RDS)
) {
  partial_errors <- readRDS(PARTIAL_ERRORS_RDS)

  if (
    is.data.frame(partial_errors) ||
      data.table::is.data.table(partial_errors)
  ) {
    comparison_errors <- data.table::as.data.table(
      partial_errors
    )
  }
}

# Keep only the expected comparison rows from a previous run.
comparison_metrics <- comparison_metrics[
  model %in% c(
    "Cox Meth+Exp",
    "Cox Meth+Exp+Clin",
    "rorS",
    "rorS+Clin"
  )
]

n_samples <- length(common_samples)
train_size <- as.integer(
  TRAINING_PROPORTION * n_samples
)

for (split_id in seq_len(N_MONTE_CARLO_SPLITS)) {
  message(
    "Comparison split ", split_id,
    "/", N_MONTE_CARLO_SPLITS
  )

  set.seed(split_id)
  train_index <- sample.int(
    n = n_samples,
    size = train_size,
    replace = FALSE
  )
  test_index <- setdiff(
    seq_len(n_samples),
    train_index
  )

  train_ids <- common_samples[train_index]
  test_ids <- common_samples[test_index]

  clin_train <- tcga_clin_matched[
    train_ids,
    ,
    drop = FALSE
  ]
  clin_test <- tcga_clin_matched[
    test_ids,
    ,
    drop = FALSE
  ]

  outcome_train <- make_survival(clin_train)

  # ----------------------------------------------------------
  # Cox-Ridge block
  # A Cox failure must not erase the rorS results.
  # ----------------------------------------------------------

  existing_cox <- comparison_metrics[
    split == split_id &
      model %in% c(
        "Cox Meth+Exp",
        "Cox Meth+Exp+Clin"
      ) &
      is.finite(c_index),
    .N
  ]

  if (existing_cox == 2L) {
    message("  Cox models already complete.")
  } else {
    cox_rows <- tryCatch(
      {
        ridge <- fit_ridge_risk(
          sample_by_feature = cox_sample_by_feature[
            train_ids,
            ,
            drop = FALSE
          ],
          outcome = outcome_train,
          split_id = split_id,
          n_total = n_samples,
          train_size = train_size
        )

        cox_train_risk <- ridge$risk

        cox_test_risk <- as.numeric(
          stats::predict(
            ridge$fit,
            newx = cox_sample_by_feature[
              test_ids,
              ,
              drop = FALSE
            ],
            s = "lambda.min",
            type = "response"
          )
        )

        cox_combined_fit <- fit_combined_model(
          molecular_risk = cox_train_risk,
          clin_train = clin_train
        )

        cox_combined_test <- predict_combined_model(
          fit = cox_combined_fit,
          molecular_risk = cox_test_risk,
          clin_test = clin_test
        )

        data.table::data.table(
          split = rep(split_id, 2L),
          cohort = rep("TCGA_test", 2L),
          model = c(
            "Cox Meth+Exp",
            "Cox Meth+Exp+Clin"
          ),
          c_index = c(
            safe_cindex(
              cox_test_risk,
              clin_test[[TIME_COL]],
              clin_test[[EVENT_COL]]
            ),
            safe_cindex(
              cox_combined_test,
              clin_test[[TIME_COL]],
              clin_test[[EVENT_COL]]
            )
          )
        )
      },
      error = function(e) {
        msg <- conditionMessage(e)

        warning(
          "Cox block failed at split ",
          split_id, ": ", msg
        )

        comparison_errors <<- comparison_errors[
          !(
            split == split_id &
              model_group == "Cox Meth+Exp"
          )
        ]

        comparison_errors <<- data.table::rbindlist(
          list(
            comparison_errors,
            data.table::data.table(
              split = split_id,
              model_group = "Cox Meth+Exp",
              error_message = msg
            )
          ),
          use.names = TRUE,
          fill = TRUE
        )

        data.table::data.table(
          split = rep(split_id, 2L),
          cohort = rep("TCGA_test", 2L),
          model = c(
            "Cox Meth+Exp",
            "Cox Meth+Exp+Clin"
          ),
          c_index = rep(NA_real_, 2L)
        )
      }
    )

    comparison_metrics <- comparison_metrics[
      !(
        split == split_id &
          model %in% c(
            "Cox Meth+Exp",
            "Cox Meth+Exp+Clin"
          )
      )
    ]

    comparison_metrics <- data.table::rbindlist(
      list(comparison_metrics, cox_rows),
      use.names = TRUE,
      fill = TRUE
    )
  }

  # ----------------------------------------------------------
  # PAM50 rorS block
  # This runs independently of Cox-Ridge.
  # ----------------------------------------------------------

  existing_rors <- comparison_metrics[
    split == split_id &
      model %in% c("rorS", "rorS+Clin") &
      is.finite(c_index),
    .N
  ]

  if (existing_rors == 2L) {
    message("  rorS models already complete.")
  } else {
    rors_rows <- tryCatch(
      {
        rors_train <- as.numeric(rors[train_ids])
        rors_test <- as.numeric(rors[test_ids])

        if (
          any(!is.finite(rors_train)) ||
            any(!is.finite(rors_test))
        ) {
          stop(
            "Non-finite rorS values after sample alignment."
          )
        }

        rors_combined_fit <- fit_combined_model(
          molecular_risk = rors_train,
          clin_train = clin_train
        )

        rors_combined_test <- predict_combined_model(
          fit = rors_combined_fit,
          molecular_risk = rors_test,
          clin_test = clin_test
        )

        data.table::data.table(
          split = rep(split_id, 2L),
          cohort = rep("TCGA_test", 2L),
          model = c("rorS", "rorS+Clin"),
          c_index = c(
            safe_cindex(
              rors_test,
              clin_test[[TIME_COL]],
              clin_test[[EVENT_COL]]
            ),
            safe_cindex(
              rors_combined_test,
              clin_test[[TIME_COL]],
              clin_test[[EVENT_COL]]
            )
          )
        )
      },
      error = function(e) {
        msg <- conditionMessage(e)

        warning(
          "rorS block failed at split ",
          split_id, ": ", msg
        )

        comparison_errors <<- comparison_errors[
          !(
            split == split_id &
              model_group == "rorS"
          )
        ]

        comparison_errors <<- data.table::rbindlist(
          list(
            comparison_errors,
            data.table::data.table(
              split = split_id,
              model_group = "rorS",
              error_message = msg
            )
          ),
          use.names = TRUE,
          fill = TRUE
        )

        data.table::data.table(
          split = rep(split_id, 2L),
          cohort = rep("TCGA_test", 2L),
          model = c("rorS", "rorS+Clin"),
          c_index = rep(NA_real_, 2L)
        )
      }
    )

    comparison_metrics <- comparison_metrics[
      !(
        split == split_id &
          model %in% c("rorS", "rorS+Clin")
      )
    ]

    comparison_metrics <- data.table::rbindlist(
      list(comparison_metrics, rors_rows),
      use.names = TRUE,
      fill = TRUE
    )
  }

  data.table::setorder(
    comparison_metrics,
    split,
    model
  )

  data.table::setorder(
    comparison_errors,
    split,
    model_group
  )

  saveRDS(
    comparison_metrics,
    PARTIAL_RDS
  )

  saveRDS(
    comparison_errors,
    PARTIAL_ERRORS_RDS
  )
}

data.table::fwrite(
  comparison_errors,
  OUTPUT_ERRORS
)

# ------------------------------------------------------------
# Combine all 9 models
# ------------------------------------------------------------

all_metrics <- data.table::rbindlist(
  list(
    full_metrics,
    exp_metrics,
    comparison_metrics
  ),
  use.names = TRUE,
  fill = TRUE
)

model_levels <- c(
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

all_metrics <- all_metrics[
  cohort == "TCGA_test" &
    model %in% model_levels
]

# Keep one row per split/model.
data.table::setorder(
  all_metrics,
  split,
  model
)

all_metrics <- unique(
  all_metrics,
  by = c("split", "model")
)

audit <- validate_model_counts(
  all_metrics,
  N_MONTE_CARLO_SPLITS
)

data.table::fwrite(
  audit,
  OUTPUT_AUDIT
)

print(audit)

if (!all(audit$complete)) {
  incomplete <- audit[complete == FALSE]

  data.table::fwrite(
    all_metrics,
    OUTPUT_METRICS
  )

  stop(
    "Not all 9 models have ",
    N_MONTE_CARLO_SPLITS,
    " valid C-index values.\nSee: ",
    OUTPUT_AUDIT,
    "\nIncomplete models:\n- ",
    paste(incomplete$model, collapse = "\n- "),
    "\n\nComparison error log: ",
    OUTPUT_ERRORS
  )
}

all_metrics[, model := factor(
  model,
  levels = model_levels,
  ordered = TRUE
)]

data.table::setorder(
  all_metrics,
  model,
  split
)

all_metrics[, model := as.character(model)]

summary_metrics <- all_metrics[
  ,
  .(
    n = .N,
    median_Cindex = stats::median(c_index, na.rm = TRUE),
    mean_Cindex = mean(c_index, na.rm = TRUE),
    Q1 = stats::quantile(c_index, 0.25, na.rm = TRUE),
    Q3 = stats::quantile(c_index, 0.75, na.rm = TRUE),
    min = min(c_index, na.rm = TRUE),
    max = max(c_index, na.rm = TRUE)
  ),
  by = model
]

summary_metrics[, order_id := match(model, model_levels)]
data.table::setorder(summary_metrics, order_id)
summary_metrics[, order_id := NULL]

data.table::fwrite(
  all_metrics,
  OUTPUT_METRICS
)

data.table::fwrite(
  summary_metrics,
  OUTPUT_SUMMARY
)

saveRDS(
  list(
    metrics = all_metrics,
    summary = summary_metrics,
    audit = audit,
    matched_samples = common_samples,
    selected_trans_genes = m2e$trans_genes,
    selected_methylation_probes = m2e$methylation_probes,
    cox_methylation_probes = selected_cox_probes,
    rors = rors
  ),
  OUTPUT_RDS
)

message("\nDONE")
message("Metrics : ", OUTPUT_METRICS)
message("Summary : ", OUTPUT_SUMMARY)
message("Audit   : ", OUTPUT_AUDIT)
message("Errors  : ", OUTPUT_ERRORS)
message("RDS     : ", OUTPUT_RDS)

print(summary_metrics)
