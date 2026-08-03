#!/usr/bin/env Rscript

# ============================================================
# STEP 09 — EXTERNAL OS DATA AUDIT FOR FIGURE 2B
#
# Purpose:
#   Audit TCGA, Terunuma, and Kao expression/clinical files before
#   running the expression-only internal/external OS models.
#
# This script DOES NOT:
#   - train any model;
#   - run ComBat/TDM;
#   - alter the input files;
#   - invent an NCA signature.
#
# Expected repository layout:
#   BRCA/
#     R/
#     results/
#     inst/extdata/csv_output/
#       TCGA_BRCA_EXP.csv
#       TCGA_BRCA_CLIN.csv
#       TERUNUMA_BRCA_EXP.csv
#       TERUNUMA_BRCA_CLIN.csv
#       KAO_BRCA_EXP.csv
#       KAO_BRCA_CLIN.csv
#
# Main outputs:
#   results/09_external_os_file_inventory.csv
#   results/09_external_os_dataset_audit.csv
#   results/09_external_os_clinical_resolution.csv
#   results/09_external_os_clinical_missingness.csv
#   results/09_external_os_gene_overlap_summary.csv
#   results/09_external_os_expression_signature_coverage.csv
#   results/09_external_os_missing_signature_genes.csv
#   results/09_external_os_expression_distribution.csv
#   results/09_external_os_readiness.csv
#   results/09_external_os_common_genes.txt
#   results/09_external_os_audit.rds
#   results/09_external_os_sessionInfo.txt
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1,
  width = 160
)

ROOT_DIR <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

DATA_DIR <- file.path(
  ROOT_DIR,
  "inst",
  "extdata",
  "csv_output"
)

RESULT_DIR <- file.path(
  ROOT_DIR,
  "results"
)

dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

required_packages <- c(
  "data.table"
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
  library(data.table)
})


# ------------------------------------------------------------
# 1. Input definitions
# ------------------------------------------------------------

cohort_files <- list(
  TCGA = list(
    expression = file.path(
      DATA_DIR,
      "TCGA_BRCA_EXP.csv"
    ),
    clinical = file.path(
      DATA_DIR,
      "TCGA_BRCA_CLIN.csv"
    )
  ),
  Terunuma = list(
    expression = file.path(
      DATA_DIR,
      "TERUNUMA_BRCA_EXP.csv"
    ),
    clinical = file.path(
      DATA_DIR,
      "TERUNUMA_BRCA_CLIN.csv"
    )
  ),
  Kao = list(
    expression = file.path(
      DATA_DIR,
      "KAO_BRCA_EXP.csv"
    ),
    clinical = file.path(
      DATA_DIR,
      "KAO_BRCA_CLIN.csv"
    )
  )
)

all_input_files <- unlist(
  lapply(
    names(cohort_files),
    function(cohort) {
      c(
        cohort_files[[cohort]]$expression,
        cohort_files[[cohort]]$clinical
      )
    }
  ),
  use.names = FALSE
)

missing_input_files <- all_input_files[
  !file.exists(all_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "Missing required input files:\n",
    paste(
      missing_input_files,
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# 2. General helpers
# ------------------------------------------------------------

canonical_name <- function(x) {
  tolower(
    gsub(
      "[^a-z0-9]",
      "",
      as.character(x)
    )
  )
}

clean_character <- function(x) {
  y <- trimws(
    as.character(x)
  )

  y[
    y %in% c(
      "",
      "NA",
      "N/A",
      "NaN",
      "--",
      ".",
      "null",
      "NULL"
    )
  ] <- NA_character_

  y
}

normalize_sample_id <- function(
    x,
    cohort
) {

  y <- clean_character(x)
  y <- chartr(
    ".",
    "-",
    y
  )

  if (identical(cohort, "TCGA")) {

    # TCGA sample-level IDs used by the project are 15 characters:
    # TCGA-XX-XXXX-01
    looks_tcga <- grepl(
      "^TCGA-",
      y,
      ignore.case = TRUE
    )

    y[looks_tcga] <- substr(
      toupper(
        y[looks_tcga]
      ),
      1,
      15
    )
  }

  y
}

resolve_column <- function(
    column_names,
    aliases
) {

  canonical_columns <- canonical_name(
    column_names
  )

  canonical_aliases <- canonical_name(
    aliases
  )

  for (alias in canonical_aliases) {
    hit <- which(
      canonical_columns == alias
    )

    if (length(hit) > 0L) {
      return(
        column_names[
          hit[[1]]
        ]
      )
    }
  }

  NA_character_
}

parse_numeric <- function(x) {
  suppressWarnings(
    as.numeric(
      clean_character(x)
    )
  )
}

parse_event <- function(x) {

  raw <- clean_character(x)
  numeric_values <- suppressWarnings(
    as.numeric(raw)
  )

  output <- rep(
    NA_real_,
    length(raw)
  )

  numeric_ok <- !is.na(
    numeric_values
  )

  output[numeric_ok] <- numeric_values[
    numeric_ok
  ]

  text <- tolower(
    raw
  )

  event_terms <- c(
    "1",
    "dead",
    "deceased",
    "death",
    "event",
    "yes",
    "true",
    "died"
  )

  censored_terms <- c(
    "0",
    "alive",
    "living",
    "censored",
    "no",
    "false"
  )

  output[
    text %in% event_terms
  ] <- 1

  output[
    text %in% censored_terms
  ] <- 0

  output
}

count_missing_character <- function(x) {
  sum(
    is.na(
      clean_character(x)
    )
  )
}

safe_fraction <- function(
    numerator,
    denominator
) {

  if (
    is.na(denominator) ||
    denominator == 0
  ) {
    return(
      NA_real_
    )
  }

  numerator / denominator
}

format_unique_examples <- function(
    x,
    max_n = 12L
) {

  values <- unique(
    clean_character(x)
  )

  values <- values[
    !is.na(values)
  ]

  paste(
    head(
      values,
      max_n
    ),
    collapse = " | "
  )
}

file_md5 <- function(path) {
  unname(
    tools::md5sum(path)
  )
}

sample_numeric_values <- function(
    data_frame,
    maximum_values = 1000000L
) {

  if (ncol(data_frame) == 0L) {
    return(
      numeric(0)
    )
  }

  total_cells <- as.double(
    nrow(data_frame)
  ) * as.double(
    ncol(data_frame)
  )

  stride <- max(
    1L,
    ceiling(
      sqrt(
        total_cells /
          maximum_values
      )
    )
  )

  row_index <- seq.int(
    1L,
    nrow(data_frame),
    by = stride
  )

  column_index <- seq.int(
    1L,
    ncol(data_frame),
    by = stride
  )

  sampled <- unlist(
    data_frame[
      row_index,
      column_index,
      drop = FALSE
    ],
    use.names = FALSE
  )

  suppressWarnings(
    as.numeric(sampled)
  )
}


# ------------------------------------------------------------
# 3. Column aliases
# ------------------------------------------------------------

column_aliases <- list(
  sample_id = c(
    "sample_id",
    "sample",
    "sampleid",
    "sample_name",
    "patient_id",
    "patient",
    "case_id",
    "id"
  ),
  expression_feature_id = c(
    "gene_id",
    "gene",
    "genesymbol",
    "gene_symbol",
    "symbol",
    "feature_id",
    "probe_id",
    "id"
  ),
  os_time = c(
    "OVERALL.SURVIVAL",
    "overall_survival",
    "overallsurvival",
    "os_time",
    "ostime",
    "survival_time",
    "survivaltime",
    "survival",
    "time"
  ),
  os_event = c(
    "overall.survival.indicator",
    "overall_survival_indicator",
    "overallsurvivalindicator",
    "os_event",
    "osevent",
    "survival_event",
    "survivalevent",
    "event",
    "status",
    "vital_status",
    "vitalstatus"
  ),
  age = c(
    "age.Dx",
    "age_dx",
    "agedx",
    "age_at_diagnosis",
    "ageatdiagnosis",
    "diagnosis_age",
    "age"
  ),
  stage = c(
    "pathologic_stage",
    "pathological_stage",
    "pathologicstage",
    "stage",
    "tumor_stage",
    "tumour_stage"
  )
)


# ------------------------------------------------------------
# 4. File inventory
# ------------------------------------------------------------

file_inventory <- do.call(
  rbind,
  lapply(
    names(cohort_files),
    function(cohort) {

      do.call(
        rbind,
        lapply(
          c(
            "expression",
            "clinical"
          ),
          function(file_type) {

            path <- cohort_files[[cohort]][[file_type]]

            info <- file.info(
              path
            )

            data.frame(
              cohort = cohort,
              file_type = file_type,
              path = normalizePath(
                path,
                winslash = "/",
                mustWork = TRUE
              ),
              size_bytes = as.numeric(
                info$size
              ),
              size_mb = as.numeric(
                info$size
              ) / 1024^2,
              modified_time = format(
                info$mtime,
                "%Y-%m-%d %H:%M:%S"
              ),
              md5 = file_md5(
                path
              ),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

write.csv(
  file_inventory,
  file.path(
    RESULT_DIR,
    "09_external_os_file_inventory.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Expression audit
# ------------------------------------------------------------

expression_objects <- list()
expression_audit_rows <- list()
expression_distribution_rows <- list()

for (cohort in names(cohort_files)) {

  expression_file <- cohort_files[[cohort]]$expression

  cat(
    "\nReading ",
    cohort,
    " expression: ",
    basename(expression_file),
    "\n",
    sep = ""
  )

  expression_data <- data.table::fread(
    expression_file,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c(
      "",
      "NA",
      "N/A",
      "NaN",
      "--"
    ),
    showProgress = TRUE
  )

  feature_column <- resolve_column(
    colnames(expression_data),
    column_aliases$expression_feature_id
  )

  if (is.na(feature_column)) {
    feature_column <- colnames(
      expression_data
    )[[1]]
  }

  feature_ids <- clean_character(
    expression_data[[feature_column]]
  )

  sample_columns <- setdiff(
    colnames(expression_data),
    feature_column
  )

  original_sample_ids <- sample_columns

  normalized_sample_ids <- normalize_sample_id(
    original_sample_ids,
    cohort
  )

  if (identical(cohort, "TCGA")) {
    primary_tumor <- grepl(
      "-01$",
      normalized_sample_ids
    )

    sample_columns <- sample_columns[
      primary_tumor
    ]

    normalized_sample_ids <- normalized_sample_ids[
      primary_tumor
    ]
  }

  molecular_values <- expression_data[
    ,
    sample_columns,
    drop = FALSE
  ]

  molecular_values[] <- lapply(
    molecular_values,
    function(column) {
      suppressWarnings(
        as.numeric(column)
      )
    }
  )

  duplicated_feature_count <- sum(
    duplicated(feature_ids),
    na.rm = TRUE
  )

  duplicated_sample_count <- sum(
    duplicated(normalized_sample_ids),
    na.rm = TRUE
  )

  missing_feature_count <- sum(
    is.na(feature_ids)
  )

  missing_cell_count <- sum(
    is.na(molecular_values)
  )

  total_cell_count <- as.double(
    nrow(molecular_values)
  ) * as.double(
    ncol(molecular_values)
  )

  sampled_values <- sample_numeric_values(
    molecular_values
  )

  sampled_values <- sampled_values[
    is.finite(sampled_values)
  ]

  if (length(sampled_values) > 0L) {
    distribution_row <- data.frame(
      cohort = cohort,
      sampled_numeric_values = length(
        sampled_values
      ),
      minimum = min(
        sampled_values
      ),
      q01 = unname(
        quantile(
          sampled_values,
          0.01
        )
      ),
      q25 = unname(
        quantile(
          sampled_values,
          0.25
        )
      ),
      median = median(
        sampled_values
      ),
      mean = mean(
        sampled_values
      ),
      q75 = unname(
        quantile(
          sampled_values,
          0.75
        )
      ),
      q99 = unname(
        quantile(
          sampled_values,
          0.99
        )
      ),
      maximum = max(
        sampled_values
      ),
      sd = sd(
        sampled_values
      ),
      fraction_negative = mean(
        sampled_values < 0
      ),
      fraction_zero = mean(
        sampled_values == 0
      ),
      stringsAsFactors = FALSE
    )
  } else {
    distribution_row <- data.frame(
      cohort = cohort,
      sampled_numeric_values = 0,
      minimum = NA_real_,
      q01 = NA_real_,
      q25 = NA_real_,
      median = NA_real_,
      mean = NA_real_,
      q75 = NA_real_,
      q99 = NA_real_,
      maximum = NA_real_,
      sd = NA_real_,
      fraction_negative = NA_real_,
      fraction_zero = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  expression_distribution_rows[[cohort]] <- distribution_row

  expression_audit_rows[[cohort]] <- data.frame(
    cohort = cohort,
    expression_rows = nrow(
      expression_data
    ),
    expression_feature_id_column = feature_column,
    expression_unique_features = length(
      unique(
        feature_ids[
          !is.na(feature_ids)
        ]
      )
    ),
    expression_missing_feature_ids = missing_feature_count,
    expression_duplicated_feature_ids = duplicated_feature_count,
    expression_samples_before_TCGA_primary_filter = length(
      original_sample_ids
    ),
    expression_samples_used = length(
      normalized_sample_ids
    ),
    expression_duplicated_sample_ids = duplicated_sample_count,
    expression_missing_cells = missing_cell_count,
    expression_total_cells = total_cell_count,
    expression_missing_fraction = safe_fraction(
      missing_cell_count,
      total_cell_count
    ),
    stringsAsFactors = FALSE
  )

  expression_objects[[cohort]] <- list(
    feature_ids = feature_ids,
    unique_feature_ids = unique(
      feature_ids[
        !is.na(feature_ids)
      ]
    ),
    sample_ids = normalized_sample_ids,
    original_sample_ids = original_sample_ids,
    duplicated_feature_ids = unique(
      feature_ids[
        duplicated(
          feature_ids
        ) &
          !is.na(feature_ids)
      ]
    ),
    duplicated_sample_ids = unique(
      normalized_sample_ids[
        duplicated(
          normalized_sample_ids
        ) &
          !is.na(normalized_sample_ids)
      ]
    )
  )

  rm(
    expression_data,
    molecular_values,
    sampled_values
  )

  invisible(
    gc()
  )
}

expression_audit <- do.call(
  rbind,
  expression_audit_rows
)

expression_distribution <- do.call(
  rbind,
  expression_distribution_rows
)

write.csv(
  expression_distribution,
  file.path(
    RESULT_DIR,
    "09_external_os_expression_distribution.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Clinical audit and sample matching
# ------------------------------------------------------------

clinical_objects <- list()
clinical_resolution_rows <- list()
clinical_missingness_rows <- list()
sample_matching_rows <- list()
clinical_audit_rows <- list()

for (cohort in names(cohort_files)) {

  clinical_file <- cohort_files[[cohort]]$clinical

  cat(
    "\nReading ",
    cohort,
    " clinical: ",
    basename(clinical_file),
    "\n",
    sep = ""
  )

  clinical_data <- data.table::fread(
    clinical_file,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c(
      "",
      "NA",
      "N/A",
      "NaN",
      "--"
    )
  )

  sample_column <- resolve_column(
    colnames(clinical_data),
    column_aliases$sample_id
  )

  if (is.na(sample_column)) {
    sample_column <- colnames(
      clinical_data
    )[[1]]
  }

  time_column <- resolve_column(
    colnames(clinical_data),
    column_aliases$os_time
  )

  event_column <- resolve_column(
    colnames(clinical_data),
    column_aliases$os_event
  )

  age_column <- resolve_column(
    colnames(clinical_data),
    column_aliases$age
  )

  stage_column <- resolve_column(
    colnames(clinical_data),
    column_aliases$stage
  )

  clinical_sample_ids <- normalize_sample_id(
    clinical_data[[sample_column]],
    cohort
  )

  clinical_duplicated_sample_ids <- unique(
    clinical_sample_ids[
      duplicated(
        clinical_sample_ids
      ) &
        !is.na(clinical_sample_ids)
    ]
  )

  clinical_resolution_rows[[cohort]] <- data.frame(
    cohort = cohort,
    sample_id_column = sample_column,
    os_time_column = time_column,
    os_event_column = event_column,
    age_column = age_column,
    stage_column = stage_column,
    all_required_columns_resolved = all(
      !is.na(
        c(
          sample_column,
          time_column,
          event_column,
          age_column,
          stage_column
        )
      )
    ),
    stringsAsFactors = FALSE
  )

  for (column_name in colnames(clinical_data)) {
    clinical_missingness_rows[[paste(cohort, column_name, sep = "::")]] <- data.frame(
      cohort = cohort,
      column = column_name,
      class = paste(
        class(
          clinical_data[[column_name]]
        ),
        collapse = "|"
      ),
      n = nrow(
        clinical_data
      ),
      n_missing = count_missing_character(
        clinical_data[[column_name]]
      ),
      missing_fraction = safe_fraction(
        count_missing_character(
          clinical_data[[column_name]]
        ),
        nrow(
          clinical_data
        )
      ),
      unique_nonmissing = length(
        unique(
          clean_character(
            clinical_data[[column_name]]
          )
        )
      ) -
        as.integer(
          any(
            is.na(
              clean_character(
                clinical_data[[column_name]]
              )
            )
          )
        ),
      examples = format_unique_examples(
        clinical_data[[column_name]]
      ),
      stringsAsFactors = FALSE
    )
  }

  time_values <- if (!is.na(time_column)) {
    parse_numeric(
      clinical_data[[time_column]]
    )
  } else {
    rep(
      NA_real_,
      nrow(clinical_data)
    )
  }

  event_values <- if (!is.na(event_column)) {
    parse_event(
      clinical_data[[event_column]]
    )
  } else {
    rep(
      NA_real_,
      nrow(clinical_data)
    )
  }

  age_values <- if (!is.na(age_column)) {
    parse_numeric(
      clinical_data[[age_column]]
    )
  } else {
    rep(
      NA_real_,
      nrow(clinical_data)
    )
  }

  stage_values <- if (!is.na(stage_column)) {
    clean_character(
      clinical_data[[stage_column]]
    )
  } else {
    rep(
      NA_character_,
      nrow(clinical_data)
    )
  }

  complete_outcome <- (
    is.finite(time_values) &
      time_values >= 0 &
      is.finite(event_values) &
      event_values %in% c(0, 1)
  )

  complete_outcome_covariates <- (
    complete_outcome &
      is.finite(age_values) &
      !is.na(stage_values)
  )

  valid_binary_event <- (
    !is.na(event_values) &
      event_values %in% c(
        0,
        1
      )
  )

  expression_sample_ids <- expression_objects[[cohort]]$sample_ids

  common_samples <- intersect(
    expression_sample_ids,
    clinical_sample_ids
  )

  expression_only_samples <- setdiff(
    expression_sample_ids,
    clinical_sample_ids
  )

  clinical_only_samples <- setdiff(
    clinical_sample_ids,
    expression_sample_ids
  )

  clinical_complete_ids <- clinical_sample_ids[
    complete_outcome_covariates
  ]

  usable_common_samples <- intersect(
    common_samples,
    clinical_complete_ids
  )

  sample_matching_rows[[cohort]] <- data.frame(
    cohort = cohort,
    expression_samples = length(
      expression_sample_ids
    ),
    clinical_rows = nrow(
      clinical_data
    ),
    clinical_unique_sample_ids = length(
      unique(
        clinical_sample_ids[
          !is.na(clinical_sample_ids)
        ]
      )
    ),
    common_expression_clinical_samples = length(
      common_samples
    ),
    expression_samples_without_clinical = length(
      expression_only_samples
    ),
    clinical_samples_without_expression = length(
      clinical_only_samples
    ),
    usable_complete_case_common_samples = length(
      usable_common_samples
    ),
    stringsAsFactors = FALSE
  )

  clinical_audit_rows[[cohort]] <- data.frame(
    cohort = cohort,
    clinical_rows = nrow(
      clinical_data
    ),
    clinical_columns = ncol(
      clinical_data
    ),
    clinical_missing_sample_ids = sum(
      is.na(
        clinical_sample_ids
      )
    ),
    clinical_duplicated_sample_ids = length(
      clinical_duplicated_sample_ids
    ),
    os_time_finite = sum(
      is.finite(
        time_values
      )
    ),
    os_time_min = if (
      any(
        is.finite(
          time_values
        )
      )
    ) {
      min(
        time_values[
          is.finite(
            time_values
          )
        ]
      )
    } else {
      NA_real_
    },
    os_time_max = if (
      any(
        is.finite(
          time_values
        )
      )
    ) {
      max(
        time_values[
          is.finite(
            time_values
          )
        ]
      )
    } else {
      NA_real_
    },
    os_event_parsed = sum(
      !is.na(
        event_values
      )
    ),
    os_event_binary_0_or_1 = sum(
      valid_binary_event
    ),
    os_events = sum(
      event_values == 1,
      na.rm = TRUE
    ),
    os_censored = sum(
      event_values == 0,
      na.rm = TRUE
    ),
    age_finite = sum(
      is.finite(
        age_values
      )
    ),
    stage_nonmissing = sum(
      !is.na(
        stage_values
      )
    ),
    complete_OS_only = sum(
      complete_outcome
    ),
    complete_OS_age_stage = sum(
      complete_outcome_covariates
    ),
    stringsAsFactors = FALSE
  )

  clinical_objects[[cohort]] <- list(
    data = clinical_data,
    sample_ids = clinical_sample_ids,
    time_column = time_column,
    event_column = event_column,
    age_column = age_column,
    stage_column = stage_column,
    time = time_values,
    event = event_values,
    age = age_values,
    stage = stage_values,
    complete_outcome = complete_outcome,
    complete_outcome_covariates = complete_outcome_covariates,
    common_samples = common_samples,
    usable_common_samples = usable_common_samples,
    expression_only_samples = expression_only_samples,
    clinical_only_samples = clinical_only_samples,
    duplicated_sample_ids = clinical_duplicated_sample_ids
  )
}

clinical_resolution <- do.call(
  rbind,
  clinical_resolution_rows
)

clinical_missingness <- do.call(
  rbind,
  clinical_missingness_rows
)

sample_matching <- do.call(
  rbind,
  sample_matching_rows
)

clinical_audit <- do.call(
  rbind,
  clinical_audit_rows
)

write.csv(
  clinical_resolution,
  file.path(
    RESULT_DIR,
    "09_external_os_clinical_resolution.csv"
  ),
  row.names = FALSE
)

write.csv(
  clinical_missingness,
  file.path(
    RESULT_DIR,
    "09_external_os_clinical_missingness.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_matching,
  file.path(
    RESULT_DIR,
    "09_external_os_sample_matching.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Combined dataset audit
# ------------------------------------------------------------

dataset_audit <- merge(
  expression_audit,
  clinical_audit,
  by = "cohort",
  all = TRUE,
  sort = FALSE
)

dataset_audit <- merge(
  dataset_audit,
  sample_matching,
  by = "cohort",
  all = TRUE,
  sort = FALSE
)

dataset_audit <- dataset_audit[
  match(
    names(cohort_files),
    dataset_audit$cohort
  ),
  ,
  drop = FALSE
]

write.csv(
  dataset_audit,
  file.path(
    RESULT_DIR,
    "09_external_os_dataset_audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Gene overlap across cohorts
# ------------------------------------------------------------

tcga_genes <- expression_objects$TCGA$unique_feature_ids
terunuma_genes <- expression_objects$Terunuma$unique_feature_ids
kao_genes <- expression_objects$Kao$unique_feature_ids

common_all_genes <- Reduce(
  intersect,
  list(
    tcga_genes,
    terunuma_genes,
    kao_genes
  )
)

gene_overlap_summary <- data.frame(
  metric = c(
    "TCGA_unique_genes",
    "Terunuma_unique_genes",
    "Kao_unique_genes",
    "TCGA_intersect_Terunuma",
    "TCGA_intersect_Kao",
    "Terunuma_intersect_Kao",
    "TCGA_intersect_Terunuma_intersect_Kao"
  ),
  value = c(
    length(
      tcga_genes
    ),
    length(
      terunuma_genes
    ),
    length(
      kao_genes
    ),
    length(
      intersect(
        tcga_genes,
        terunuma_genes
      )
    ),
    length(
      intersect(
        tcga_genes,
        kao_genes
      )
    ),
    length(
      intersect(
        terunuma_genes,
        kao_genes
      )
    ),
    length(
      common_all_genes
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  gene_overlap_summary,
  file.path(
    RESULT_DIR,
    "09_external_os_gene_overlap_summary.csv"
  ),
  row.names = FALSE
)

writeLines(
  sort(
    common_all_genes
  ),
  file.path(
    RESULT_DIR,
    "09_external_os_common_genes.txt"
  )
)


# ------------------------------------------------------------
# 9. Load this replication's expression-only M2EFM signature
# ------------------------------------------------------------

signature_candidates <- c(
  file.path(
    RESULT_DIR,
    "02_expression_only_m2eGenes.txt"
  ),
  file.path(
    ROOT_DIR,
    "02_expression_only_m2eGenes.txt"
  )
)

signature_file <- signature_candidates[
  file.exists(
    signature_candidates
  )
]

signature_source <- NA_character_
expression_signature <- character(0)

if (length(signature_file) > 0L) {

  signature_source <- normalizePath(
    signature_file[[1]],
    winslash = "/",
    mustWork = TRUE
  )

  expression_signature <- unique(
    clean_character(
      readLines(
        signature_file[[1]],
        warn = FALSE
      )
    )
  )

  expression_signature <- expression_signature[
    !is.na(
      expression_signature
    )
  ]

} else {

  m2eqtl_file <- file.path(
    RESULT_DIR,
    "01_TCGA_OS_m2eQTL.rds"
  )

  m2efm_source <- file.path(
    ROOT_DIR,
    "R",
    "m2efm.R"
  )

  if (
    file.exists(
      m2eqtl_file
    ) &&
      file.exists(
        m2efm_source
      )
  ) {

    source(
      m2efm_source,
      local = .GlobalEnv
    )

    m2eqtl_object <- readRDS(
      m2eqtl_file
    )

    expression_signature <- unique(
      as.character(
        get_genes(
          m2eqtl_object,
          gene_type = "trans",
          integrate_data = FALSE
        )
      )
    )

    signature_source <- paste0(
      normalizePath(
        m2eqtl_file,
        winslash = "/",
        mustWork = TRUE
      ),
      " via get_genes(..., integrate_data=FALSE)"
    )
  }
}

signature_coverage_rows <- list()
missing_signature_rows <- list()

if (length(expression_signature) > 0L) {

  for (cohort in names(expression_objects)) {

    cohort_genes <- expression_objects[[cohort]]$unique_feature_ids

    present <- expression_signature %in%
      cohort_genes

    signature_coverage_rows[[cohort]] <- data.frame(
      cohort = cohort,
      signature_source = signature_source,
      signature_genes = length(
        expression_signature
      ),
      genes_present = sum(
        present
      ),
      genes_missing = sum(
        !present
      ),
      coverage_fraction = mean(
        present
      ),
      stringsAsFactors = FALSE
    )

    missing_genes_for_cohort <- expression_signature[
      !present
    ]

    missing_signature_rows[[cohort]] <- data.frame(
      cohort = rep(
        cohort,
        length(missing_genes_for_cohort)
      ),
      gene = missing_genes_for_cohort,
      stringsAsFactors = FALSE
    )
  }

  all_present <- expression_signature %in%
    common_all_genes

  signature_coverage_rows[["Common_all"]] <- data.frame(
    cohort = "Common_all",
    signature_source = signature_source,
    signature_genes = length(
      expression_signature
    ),
    genes_present = sum(
      all_present
    ),
    genes_missing = sum(
      !all_present
    ),
    coverage_fraction = mean(
      all_present
    ),
    stringsAsFactors = FALSE
  )

  missing_genes_common_all <- expression_signature[
    !all_present
  ]

  missing_signature_rows[["Common_all"]] <- data.frame(
    cohort = rep(
      "Common_all",
      length(missing_genes_common_all)
    ),
    gene = missing_genes_common_all,
    stringsAsFactors = FALSE
  )

  signature_coverage <- do.call(
    rbind,
    signature_coverage_rows
  )

  missing_signature_genes <- do.call(
    rbind,
    missing_signature_rows
  )

} else {

  signature_coverage <- data.frame(
    cohort = c(
      names(expression_objects),
      "Common_all"
    ),
    signature_source = NA_character_,
    signature_genes = 0,
    genes_present = 0,
    genes_missing = NA_integer_,
    coverage_fraction = NA_real_,
    stringsAsFactors = FALSE
  )

  missing_signature_genes <- data.frame(
    cohort = NA_character_,
    gene = NA_character_,
    stringsAsFactors = FALSE
  )
}

write.csv(
  signature_coverage,
  file.path(
    RESULT_DIR,
    "09_external_os_expression_signature_coverage.csv"
  ),
  row.names = FALSE
)

write.csv(
  missing_signature_genes,
  file.path(
    RESULT_DIR,
    "09_external_os_missing_signature_genes.csv"
  ),
  row.names = FALSE
)

if (length(expression_signature) > 0L) {
  writeLines(
    expression_signature,
    file.path(
      RESULT_DIR,
      "09_external_os_expression_signature.txt"
    )
  )
}


# ------------------------------------------------------------
# 10. Search for an explicitly supplied NCA signature
# ------------------------------------------------------------

all_project_files <- list.files(
  ROOT_DIR,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

nca_candidate_files <- all_project_files[
  basename(all_project_files) !=
    "09_external_os_nca_signature_audit.csv" &
  grepl(
    "nca",
    basename(
      all_project_files
    ),
    ignore.case = TRUE
  ) &
    grepl(
      "gene|signature|list",
      basename(
        all_project_files
      ),
      ignore.case = TRUE
    ) &
    grepl(
      "\\.(txt|csv|tsv)$",
      basename(
        all_project_files
      ),
      ignore.case = TRUE
    )
]

nca_audit <- data.frame(
  nca_signature_found = length(
    nca_candidate_files
  ) > 0L,
  candidate_count = length(
    nca_candidate_files
  ),
  candidate_files = paste(
    normalizePath(
      nca_candidate_files,
      winslash = "/",
      mustWork = FALSE
    ),
    collapse = " | "
  ),
  note = if (
    length(
      nca_candidate_files
    ) > 0L
  ) {
    paste(
      "Candidate file(s) found. Content must be manually verified",
      "before calling the models NCA."
    )
  } else {
    paste(
      "No explicit NCA gene-signature file found.",
      "Do not reconstruct NCA models by inventing a list."
    )
  },
  stringsAsFactors = FALSE
)

write.csv(
  nca_audit,
  file.path(
    RESULT_DIR,
    "09_external_os_nca_signature_audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Readiness checks
# ------------------------------------------------------------

resolution_by_cohort <- setNames(
  clinical_resolution$all_required_columns_resolved,
  clinical_resolution$cohort
)

usable_samples_by_cohort <- setNames(
  sample_matching$usable_complete_case_common_samples,
  sample_matching$cohort
)

signature_common_row <- signature_coverage[
  signature_coverage$cohort == "Common_all",
  ,
  drop = FALSE
]

signature_ready <- (
  nrow(
    signature_common_row
  ) == 1L &&
    signature_common_row$signature_genes > 0L &&
    signature_common_row$genes_missing == 0L
)

readiness <- data.frame(
  check = c(
    "all_six_required_input_files_exist",
    "TCGA_required_clinical_columns_resolved",
    "Terunuma_required_clinical_columns_resolved",
    "Kao_required_clinical_columns_resolved",
    "TCGA_has_usable_complete_case_samples",
    "Terunuma_has_usable_complete_case_samples",
    "Kao_has_usable_complete_case_samples",
    "three_way_common_gene_set_nonempty",
    "replication_expression_signature_loaded",
    "replication_expression_signature_present_in_all_three_cohorts",
    "explicit_NCA_signature_candidate_found",
    "batch_harmonization_status_confirmed"
  ),
  status = c(
    all(
      file.exists(
        all_input_files
      )
    ),
    isTRUE(
      resolution_by_cohort[["TCGA"]]
    ),
    isTRUE(
      resolution_by_cohort[["Terunuma"]]
    ),
    isTRUE(
      resolution_by_cohort[["Kao"]]
    ),
    usable_samples_by_cohort[["TCGA"]] > 0,
    usable_samples_by_cohort[["Terunuma"]] > 0,
    usable_samples_by_cohort[["Kao"]] > 0,
    length(
      common_all_genes
    ) > 0L,
    length(
      expression_signature
    ) > 0L,
    signature_ready,
    nca_audit$nca_signature_found,
    FALSE
  ),
  severity = c(
    "required",
    "required",
    "required",
    "required",
    "required",
    "required",
    "required",
    "required",
    "required",
    "required",
    "optional_for_7_models_required_for_9_models",
    "manual_review"
  ),
  note = c(
    "All expected processed CSV files must be present.",
    "Requires sample ID, OS time, OS event, age, and stage.",
    "Requires sample ID, OS time, OS event, age, and stage.",
    "Requires sample ID, OS time, OS event, age, and stage.",
    "Full expression-only TCGA cohort after complete-case filtering.",
    "External validation cohort after complete-case filtering.",
    "External validation cohort after complete-case filtering.",
    "Required for common-gene Cox expression baseline and harmonization.",
    paste0(
      "Loaded signature size: ",
      length(
        expression_signature
      )
    ),
    paste0(
      "Missing from three-way common genes: ",
      if (
        nrow(
          signature_common_row
        ) == 1L
      ) {
        signature_common_row$genes_missing
      } else {
        NA
      }
    ),
    nca_audit$note,
    paste(
      "This audit reports value distributions but cannot prove",
      "whether TDM/ComBat has already been applied.",
      "Decide before Step 10 to avoid double correction."
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  readiness,
  file.path(
    RESULT_DIR,
    "09_external_os_readiness.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Save an R object for the next step
# ------------------------------------------------------------

audit_object <- list(
  root_dir = ROOT_DIR,
  data_dir = DATA_DIR,
  result_dir = RESULT_DIR,
  file_inventory = file_inventory,
  expression_audit = expression_audit,
  expression_distribution = expression_distribution,
  clinical_resolution = clinical_resolution,
  clinical_missingness = clinical_missingness,
  clinical_audit = clinical_audit,
  sample_matching = sample_matching,
  dataset_audit = dataset_audit,
  gene_overlap_summary = gene_overlap_summary,
  common_all_genes = common_all_genes,
  expression_signature_source = signature_source,
  expression_signature = expression_signature,
  signature_coverage = signature_coverage,
  missing_signature_genes = missing_signature_genes,
  nca_audit = nca_audit,
  readiness = readiness,
  clinical_objects = lapply(
    clinical_objects,
    function(object) {
      object[
        setdiff(
          names(object),
          "data"
        )
      ]
    }
  ),
  expression_objects = expression_objects,
  created_at = Sys.time()
)

saveRDS(
  audit_object,
  file.path(
    RESULT_DIR,
    "09_external_os_audit.rds"
  ),
  compress = "xz"
)


# ------------------------------------------------------------
# 13. Reproducibility report
# ------------------------------------------------------------

capture.output(
  {
    cat("STEP 09 — EXTERNAL OS DATA AUDIT\n")
    cat("Created: ", format(Sys.time()), "\n\n", sep = "")

    cat("FILE INVENTORY\n")
    print(file_inventory)
    cat("\nDATASET AUDIT\n")
    print(dataset_audit)
    cat("\nCLINICAL COLUMN RESOLUTION\n")
    print(clinical_resolution)
    cat("\nSAMPLE MATCHING\n")
    print(sample_matching)
    cat("\nGENE OVERLAP\n")
    print(gene_overlap_summary)
    cat("\nEXPRESSION-ONLY SIGNATURE COVERAGE\n")
    print(signature_coverage)
    cat("\nNCA SIGNATURE AUDIT\n")
    print(nca_audit)
    cat("\nREADINESS\n")
    print(readiness)
    cat("\nEXPRESSION DISTRIBUTIONS\n")
    print(expression_distribution)
    cat("\nSESSION INFO\n")
    sessionInfo()
  },
  file = file.path(
    RESULT_DIR,
    "09_external_os_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 14. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 09 EXTERNAL OS DATA AUDIT COMPLETED\n")
cat("============================================================\n")

cat("\nUSABLE COMPLETE-CASE EXPRESSION + CLINICAL SAMPLES\n")
print(
  sample_matching[
    ,
    c(
      "cohort",
      "common_expression_clinical_samples",
      "usable_complete_case_common_samples"
    )
  ],
  row.names = FALSE
)

cat("\nTHREE-WAY GENE OVERLAP\n")
print(
  gene_overlap_summary,
  row.names = FALSE
)

cat("\nEXPRESSION-ONLY SIGNATURE COVERAGE\n")
print(
  signature_coverage,
  row.names = FALSE
)

cat("\nCLINICAL COLUMN RESOLUTION\n")
print(
  clinical_resolution,
  row.names = FALSE
)

cat("\nEXPRESSION VALUE DISTRIBUTIONS\n")
print(
  expression_distribution,
  row.names = FALSE
)

cat("\nREADINESS CHECKS\n")
print(
  readiness,
  row.names = FALSE
)

cat(
  "\nImportant: batch_harmonization_status_confirmed is intentionally FALSE.\n"
)

cat(
  "Review the distribution table and preprocessing provenance before Step 10.\n"
)

cat(
  "\nResults saved to: ",
  RESULT_DIR,
  "\n",
  sep = ""
)

cat("\nFinished successfully.\n")
