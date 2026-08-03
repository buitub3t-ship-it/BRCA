#!/usr/bin/env Rscript

# ============================================================
# STEP 09B — PREPARE TERUNUMA GENE-LEVEL EXPRESSION AND
#            RE-AUDIT FIGURE 2B INPUTS
#
# This script:
#   1. Maps Terunuma GPL6244 transcript-cluster IDs to gene symbols
#      with hugene10sttranscriptcluster.db.
#   2. Removes transcript clusters mapping ambiguously to >1 symbol.
#   3. Collapses multiple transcript clusters per gene by sample-wise median.
#   4. Locks OS to:
#        OVERALL.SURVIVAL
#        overall.survival.indicator
#   5. Aligns expression and clinical samples.
#   6. Reduces stage to Stage I/II/III/IV.
#   7. Re-audits common genes and the 119-gene expression-only signature.
#
# It DOES NOT overwrite the original Terunuma CSV.
# It DOES NOT run ComBat or train models.
#
# Main outputs:
#   results/09B_TERUNUMA_BRCA_EXP_GENE_MEDIAN.csv
#   results/09B_terunuma_gene_level_expression.rds
#   results/09B_external_os_preprocessed.rds
#   results/09B_terunuma_mapping_audit.csv
#   results/09B_terunuma_probe_to_gene_mapping.csv
#   results/09B_external_os_sample_audit.csv
#   results/09B_external_os_gene_overlap.csv
#   results/09B_external_os_signature_coverage.csv
#   results/09B_external_os_missing_signature_genes.csv
#   results/09B_external_os_expression_distribution.csv
#   results/09B_external_os_readiness.csv
#   results/09B_external_os_common_genes.txt
#   results/09B_external_os_sessionInfo.txt
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
  "data.table",
  "AnnotationDbi",
  "hugene10sttranscriptcluster.db",
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
  library(data.table)
  library(AnnotationDbi)
  library(hugene10sttranscriptcluster.db)
  library(matrixStats)
})


# ------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------

FILES <- c(
  TCGA_EXP = file.path(
    DATA_DIR,
    "TCGA_BRCA_EXP.csv"
  ),
  TCGA_CLIN = file.path(
    DATA_DIR,
    "TCGA_BRCA_CLIN.csv"
  ),
  TERU_EXP = file.path(
    DATA_DIR,
    "TERUNUMA_BRCA_EXP.csv"
  ),
  TERU_CLIN = file.path(
    DATA_DIR,
    "TERUNUMA_BRCA_CLIN.csv"
  ),
  KAO_EXP = file.path(
    DATA_DIR,
    "KAO_BRCA_EXP.csv"
  ),
  KAO_CLIN = file.path(
    DATA_DIR,
    "KAO_BRCA_CLIN.csv"
  ),
  SIGNATURE = file.path(
    RESULT_DIR,
    "02_expression_only_m2eGenes.txt"
  )
)

missing_files <- FILES[
  !file.exists(FILES)
]

if (length(missing_files) > 0L) {
  stop(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

TIME_COL <- "OVERALL.SURVIVAL"
EVENT_COL <- "overall.survival.indicator"
AGE_COL <- "age.Dx"
STAGE_COL <- "pathologic_stage"

REQUIRED_CLINICAL <- c(
  TIME_COL,
  EVENT_COL,
  AGE_COL,
  STAGE_COL
)


# ------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------

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
      "[Not Available]",
      "[Not Applicable]",
      "[Unknown]"
    )
  ] <- NA_character_

  y
}

force_numeric <- function(x) {
  suppressWarnings(
    as.numeric(
      clean_character(x)
    )
  )
}

normalize_sample_id <- function(
    x,
    cohort
) {

  y <- chartr(
    ".",
    "-",
    clean_character(x)
  )

  if (identical(cohort, "TCGA")) {
    is_tcga <- grepl(
      "^TCGA-",
      y,
      ignore.case = TRUE
    )

    y[is_tcga] <- substr(
      toupper(
        y[is_tcga]
      ),
      1,
      15
    )
  }

  y
}

reduce_stage_local <- function(x) {
  y <- toupper(
    trimws(
      clean_character(x)
    )
  )

  y <- gsub(
    "^STAGE[ _-]*",
    "",
    y
  )

  output <- rep(
    NA_character_,
    length(y)
  )

  output[
    grepl(
      "^I($|A$|B$|C$|[0-9])",
      y
    )
  ] <- "Stage I"

  output[
    grepl(
      "^II($|A$|B$|C$|[0-9])",
      y
    )
  ] <- "Stage II"

  output[
    grepl(
      "^III($|A$|B$|C$|[0-9])",
      y
    )
  ] <- "Stage III"

  output[
    grepl(
      "^IV($|A$|B$|C$|[0-9])",
      y
    )
  ] <- "Stage IV"

  # The order above can let "III" or "IV" match "^I".
  # Re-assert longer Roman numerals last.
  output[
    grepl(
      "^II($|A$|B$|C$|[0-9])",
      y
    )
  ] <- "Stage II"

  output[
    grepl(
      "^III($|A$|B$|C$|[0-9])",
      y
    )
  ] <- "Stage III"

  output[
    grepl(
      "^IV($|A$|B$|C$|[0-9])",
      y
    )
  ] <- "Stage IV"

  factor(
    output,
    levels = c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  )
}

read_profile_matrix <- function(
    path,
    cohort
) {

  cat(
    "Reading ",
    cohort,
    " expression: ",
    basename(path),
    "\n",
    sep = ""
  )

  x <- data.table::fread(
    path,
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

  if (ncol(x) < 2L) {
    stop(
      "Expression file has fewer than two columns: ",
      path
    )
  }

  feature_ids <- clean_character(
    x[[1L]]
  )

  if (
    anyNA(feature_ids) ||
    any(feature_ids == "")
  ) {
    stop(
      "Missing feature IDs in ",
      path
    )
  }

  if (anyDuplicated(feature_ids)) {
    stop(
      "Duplicated feature IDs in ",
      path
    )
  }

  x <- x[
    -1L
  ]

  matrix_values <- data.matrix(
    x
  )

  rownames(
    matrix_values
  ) <- feature_ids

  colnames(
    matrix_values
  ) <- normalize_sample_id(
    colnames(matrix_values),
    cohort
  )

  if (identical(cohort, "TCGA")) {
    primary <- grepl(
      "-01$",
      colnames(matrix_values)
    )

    matrix_values <- matrix_values[
      ,
      primary,
      drop = FALSE
    ]
  }

  if (anyDuplicated(colnames(matrix_values))) {
    stop(
      "Duplicated sample IDs after normalization in ",
      cohort
    )
  }

  if (anyNA(matrix_values)) {
    stop(
      "Expression matrix contains NA/non-numeric values in ",
      cohort
    )
  }

  matrix_values
}

read_clinical <- function(
    path,
    cohort
) {

  cat(
    "Reading ",
    cohort,
    " clinical: ",
    basename(path),
    "\n",
    sep = ""
  )

  x <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c(
      "",
      "NA",
      "N/A",
      "NaN",
      "--",
      "[Not Available]",
      "[Not Applicable]",
      "[Unknown]"
    )
  )

  if (!"sample_id" %in% colnames(x)) {
    stop(
      cohort,
      " clinical file lacks sample_id."
    )
  }

  absent <- setdiff(
    REQUIRED_CLINICAL,
    colnames(x)
  )

  if (length(absent) > 0L) {
    stop(
      cohort,
      " clinical file lacks: ",
      paste(
        absent,
        collapse = ", "
      )
    )
  }

  sample_ids <- normalize_sample_id(
    x$sample_id,
    cohort
  )

  if (
    anyNA(sample_ids) ||
    any(sample_ids == "")
  ) {
    stop(
      cohort,
      " clinical contains missing sample IDs."
    )
  }

  if (anyDuplicated(sample_ids)) {
    stop(
      cohort,
      " clinical contains duplicated sample IDs."
    )
  }

  rownames(x) <- sample_ids
  x$sample_id <- NULL

  x[[TIME_COL]] <- force_numeric(
    x[[TIME_COL]]
  )

  x[[EVENT_COL]] <- force_numeric(
    x[[EVENT_COL]]
  )

  x[[AGE_COL]] <- force_numeric(
    x[[AGE_COL]]
  )

  x[[STAGE_COL]] <- reduce_stage_local(
    x[[STAGE_COL]]
  )

  observed_events <- sort(
    unique(
      stats::na.omit(
        x[[EVENT_COL]]
      )
    )
  )

  if (!all(observed_events %in% c(0, 1))) {
    stop(
      cohort,
      " OS event must be coded 0/1. Observed: ",
      paste(
        observed_events,
        collapse = ", "
      )
    )
  }

  x
}

align_expression_clinical <- function(
    expression,
    clinical,
    cohort
) {

  common <- colnames(expression)[
    colnames(expression) %in%
      rownames(clinical)
  ]

  if (length(common) == 0L) {
    stop(
      "No matching expression/clinical samples for ",
      cohort
    )
  }

  expression <- expression[
    ,
    common,
    drop = FALSE
  ]

  clinical <- clinical[
    common,
    ,
    drop = FALSE
  ]

  stopifnot(
    identical(
      colnames(expression),
      rownames(clinical)
    )
  )

  complete <- complete.cases(
    clinical[
      ,
      REQUIRED_CLINICAL,
      drop = FALSE
    ]
  )

  complete <- complete &
    is.finite(
      clinical[[TIME_COL]]
    ) &
    is.finite(
      clinical[[EVENT_COL]]
    ) &
    is.finite(
      clinical[[AGE_COL]]
    ) &
    clinical[[TIME_COL]] >= 0

  expression_complete <- expression[
    ,
    complete,
    drop = FALSE
  ]

  clinical_complete <- clinical[
    complete,
    ,
    drop = FALSE
  ]

  clinical_complete[[STAGE_COL]] <- droplevels(
    clinical_complete[[STAGE_COL]]
  )

  list(
    expression_all_matched = expression,
    clinical_all_matched = clinical,
    expression_complete = expression_complete,
    clinical_complete = clinical_complete,
    common_samples = common,
    complete_mask = complete
  )
}

sample_distribution <- function(
    matrix_values,
    cohort,
    maximum_values = 1000000L
) {

  total_cells <- as.double(
    length(matrix_values)
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

  rows <- seq.int(
    1L,
    nrow(matrix_values),
    by = stride
  )

  columns <- seq.int(
    1L,
    ncol(matrix_values),
    by = stride
  )

  values <- as.numeric(
    matrix_values[
      rows,
      columns,
      drop = FALSE
    ]
  )

  values <- values[
    is.finite(values)
  ]

  data.frame(
    cohort = cohort,
    sampled_values = length(values),
    minimum = min(values),
    q01 = unname(
      quantile(
        values,
        0.01
      )
    ),
    q25 = unname(
      quantile(
        values,
        0.25
      )
    ),
    median = median(values),
    mean = mean(values),
    q75 = unname(
      quantile(
        values,
        0.75
      )
    ),
    q99 = unname(
      quantile(
        values,
        0.99
      )
    ),
    maximum = max(values),
    sd = sd(values),
    fraction_negative = mean(
      values < 0
    ),
    fraction_zero = mean(
      values == 0
    ),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 3. Load raw matrices and clinical tables
# ------------------------------------------------------------

tcga_exp <- read_profile_matrix(
  FILES[["TCGA_EXP"]],
  "TCGA"
)

teru_probe_exp <- read_profile_matrix(
  FILES[["TERU_EXP"]],
  "Terunuma"
)

kao_exp <- read_profile_matrix(
  FILES[["KAO_EXP"]],
  "Kao"
)

tcga_clin <- read_clinical(
  FILES[["TCGA_CLIN"]],
  "TCGA"
)

teru_clin <- read_clinical(
  FILES[["TERU_CLIN"]],
  "Terunuma"
)

kao_clin <- read_clinical(
  FILES[["KAO_CLIN"]],
  "Kao"
)


# ------------------------------------------------------------
# 4. Map Terunuma transcript clusters to gene symbols
# ------------------------------------------------------------

cat(
  "\nMapping Terunuma transcript clusters to gene symbols...\n"
)

probe_ids <- rownames(
  teru_probe_exp
)

annotation <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = probe_ids,
  keytype = "PROBEID",
  columns = c(
    "SYMBOL",
    "ENTREZID"
  )
)

annotation$PROBEID <- clean_character(
  annotation$PROBEID
)

annotation$SYMBOL <- clean_character(
  annotation$SYMBOL
)

annotation$ENTREZID <- clean_character(
  annotation$ENTREZID
)

annotation <- annotation[
  !is.na(annotation$PROBEID) &
    annotation$PROBEID %in% probe_ids &
    !is.na(annotation$SYMBOL),
  ,
  drop = FALSE
]

annotation <- unique(
  annotation[
    ,
    c(
      "PROBEID",
      "SYMBOL",
      "ENTREZID"
    )
  ]
)

# Count unique symbols per probe ID.
symbol_count <- tapply(
  annotation$SYMBOL,
  annotation$PROBEID,
  function(values) {
    length(
      unique(values)
    )
  }
)

unambiguous_probe_ids <- names(
  symbol_count[
    symbol_count == 1L
  ]
)

ambiguous_probe_ids <- names(
  symbol_count[
    symbol_count > 1L
  ]
)

mapping <- annotation[
  annotation$PROBEID %in%
    unambiguous_probe_ids,
  ,
  drop = FALSE
]

# One unambiguous probe may still have duplicate rows due repeated Entrez
# records. Keep one probe-symbol pair.
mapping <- unique(
  mapping[
    ,
    c(
      "PROBEID",
      "SYMBOL"
    )
  ]
)

if (anyDuplicated(mapping$PROBEID)) {
  stop(
    "Internal error: unambiguous mapping still contains duplicated probe IDs."
  )
}

mapped_probe_ids <- mapping$PROBEID

probe_index <- match(
  mapped_probe_ids,
  rownames(teru_probe_exp)
)

if (anyNA(probe_index)) {
  stop(
    "Some mapped probe IDs were not found in the Terunuma matrix."
  )
}

mapped_matrix <- teru_probe_exp[
  probe_index,
  ,
  drop = FALSE
]

rownames(mapped_matrix) <- mapping$SYMBOL

gene_groups <- split(
  seq_len(
    nrow(mapped_matrix)
  ),
  rownames(mapped_matrix)
)

cat(
  "Collapsing ",
  nrow(mapped_matrix),
  " unambiguous transcript clusters into ",
  length(gene_groups),
  " gene symbols by sample-wise median...\n",
  sep = ""
)

teru_gene_exp <- t(
  vapply(
    gene_groups,
    function(index) {
      matrixStats::colMedians(
        mapped_matrix[
          index,
          ,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    },
    numeric(
      ncol(mapped_matrix)
    )
  )
)

colnames(teru_gene_exp) <- colnames(
  mapped_matrix
)

if (
  anyNA(teru_gene_exp) ||
  any(!is.finite(teru_gene_exp))
) {
  stop(
    "Terunuma gene-level matrix contains missing/non-finite values."
  )
}

if (anyDuplicated(rownames(teru_gene_exp))) {
  stop(
    "Terunuma gene-level matrix still contains duplicated gene symbols."
  )
}

mapping_audit <- data.frame(
  metric = c(
    "input_transcript_clusters",
    "IDs_present_in_annotation_DB",
    "probe_IDs_with_any_gene_symbol",
    "ambiguous_probe_IDs_removed",
    "unambiguous_probe_IDs_retained",
    "output_unique_gene_symbols",
    "unmapped_or_symbol_missing_probe_IDs",
    "mapping_package",
    "mapping_package_version"
  ),
  value = c(
    length(probe_ids),
    sum(
      probe_ids %in%
        AnnotationDbi::keys(
          hugene10sttranscriptcluster.db,
          keytype = "PROBEID"
        )
    ),
    length(
      unique(
        annotation$PROBEID
      )
    ),
    length(
      ambiguous_probe_ids
    ),
    length(
      mapped_probe_ids
    ),
    nrow(
      teru_gene_exp
    ),
    length(
      setdiff(
        probe_ids,
        unique(
          annotation$PROBEID
        )
      )
    ),
    "hugene10sttranscriptcluster.db",
    as.character(
      packageVersion(
        "hugene10sttranscriptcluster.db"
      )
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  mapping_audit,
  file.path(
    RESULT_DIR,
    "09B_terunuma_mapping_audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  mapping,
  file.path(
    RESULT_DIR,
    "09B_terunuma_probe_to_gene_mapping.csv"
  ),
  row.names = FALSE
)

writeLines(
  ambiguous_probe_ids,
  file.path(
    RESULT_DIR,
    "09B_terunuma_ambiguous_probe_ids_removed.txt"
  )
)

writeLines(
  setdiff(
    probe_ids,
    unique(
      annotation$PROBEID
    )
  ),
  file.path(
    RESULT_DIR,
    "09B_terunuma_unmapped_probe_ids.txt"
  )
)

teru_gene_csv <- data.frame(
  gene_id = rownames(
    teru_gene_exp
  ),
  teru_gene_exp,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  teru_gene_csv,
  file.path(
    RESULT_DIR,
    "09B_TERUNUMA_BRCA_EXP_GENE_MEDIAN.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  teru_gene_exp,
  file.path(
    RESULT_DIR,
    "09B_terunuma_gene_level_expression.rds"
  ),
  compress = "xz"
)


# ------------------------------------------------------------
# 5. Align expression and clinical data with exact OS columns
# ------------------------------------------------------------

tcga_ready <- align_expression_clinical(
  tcga_exp,
  tcga_clin,
  "TCGA"
)

teru_ready <- align_expression_clinical(
  teru_gene_exp,
  teru_clin,
  "Terunuma"
)

kao_ready <- align_expression_clinical(
  kao_exp,
  kao_clin,
  "Kao"
)

sample_audit <- do.call(
  rbind,
  list(
    data.frame(
      cohort = "TCGA",
      expression_samples = ncol(tcga_exp),
      clinical_samples = nrow(tcga_clin),
      common_expression_clinical_samples = length(
        tcga_ready$common_samples
      ),
      samples_removed_missing_OS_age_stage = sum(
        !tcga_ready$complete_mask
      ),
      usable_complete_case_samples = ncol(
        tcga_ready$expression_complete
      ),
      OS_events = sum(
        tcga_ready$clinical_complete[[EVENT_COL]] == 1
      ),
      OS_censored = sum(
        tcga_ready$clinical_complete[[EVENT_COL]] == 0
      ),
      time_column = TIME_COL,
      event_column = EVENT_COL
    ),
    data.frame(
      cohort = "Terunuma",
      expression_samples = ncol(teru_gene_exp),
      clinical_samples = nrow(teru_clin),
      common_expression_clinical_samples = length(
        teru_ready$common_samples
      ),
      samples_removed_missing_OS_age_stage = sum(
        !teru_ready$complete_mask
      ),
      usable_complete_case_samples = ncol(
        teru_ready$expression_complete
      ),
      OS_events = sum(
        teru_ready$clinical_complete[[EVENT_COL]] == 1
      ),
      OS_censored = sum(
        teru_ready$clinical_complete[[EVENT_COL]] == 0
      ),
      time_column = TIME_COL,
      event_column = EVENT_COL
    ),
    data.frame(
      cohort = "Kao",
      expression_samples = ncol(kao_exp),
      clinical_samples = nrow(kao_clin),
      common_expression_clinical_samples = length(
        kao_ready$common_samples
      ),
      samples_removed_missing_OS_age_stage = sum(
        !kao_ready$complete_mask
      ),
      usable_complete_case_samples = ncol(
        kao_ready$expression_complete
      ),
      OS_events = sum(
        kao_ready$clinical_complete[[EVENT_COL]] == 1
      ),
      OS_censored = sum(
        kao_ready$clinical_complete[[EVENT_COL]] == 0
      ),
      time_column = TIME_COL,
      event_column = EVENT_COL
    )
  )
)

write.csv(
  sample_audit,
  file.path(
    RESULT_DIR,
    "09B_external_os_sample_audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Gene overlap and expression-only signature coverage
# ------------------------------------------------------------

tcga_genes <- rownames(
  tcga_ready$expression_complete
)

teru_genes <- rownames(
  teru_ready$expression_complete
)

kao_genes <- rownames(
  kao_ready$expression_complete
)

common_all_genes <- Reduce(
  intersect,
  list(
    tcga_genes,
    teru_genes,
    kao_genes
  )
)

# Preserve TCGA gene order.
common_all_genes <- tcga_genes[
  tcga_genes %in%
    common_all_genes
]

gene_overlap <- data.frame(
  metric = c(
    "TCGA_unique_genes",
    "Terunuma_gene_level_unique_genes",
    "Kao_unique_genes",
    "TCGA_intersect_Terunuma",
    "TCGA_intersect_Kao",
    "Terunuma_intersect_Kao",
    "three_way_common_genes"
  ),
  value = c(
    length(tcga_genes),
    length(teru_genes),
    length(kao_genes),
    length(
      intersect(
        tcga_genes,
        teru_genes
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
        teru_genes,
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
  gene_overlap,
  file.path(
    RESULT_DIR,
    "09B_external_os_gene_overlap.csv"
  ),
  row.names = FALSE
)

writeLines(
  common_all_genes,
  file.path(
    RESULT_DIR,
    "09B_external_os_common_genes.txt"
  )
)

expression_signature <- unique(
  clean_character(
    readLines(
      FILES[["SIGNATURE"]],
      warn = FALSE
    )
  )
)

expression_signature <- expression_signature[
  !is.na(expression_signature)
]

signature_sets <- list(
  TCGA = tcga_genes,
  Terunuma = teru_genes,
  Kao = kao_genes,
  Common_all = common_all_genes
)

signature_coverage <- do.call(
  rbind,
  lapply(
    names(signature_sets),
    function(cohort) {
      present <- expression_signature %in%
        signature_sets[[cohort]]

      data.frame(
        cohort = cohort,
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
    }
  )
)

missing_signature_list <- lapply(
  names(signature_sets),
  function(cohort) {
    missing_genes <- setdiff(
      expression_signature,
      signature_sets[[cohort]]
    )

    if (length(missing_genes) == 0L) {
      return(
        data.frame(
          cohort = character(0),
          gene = character(0),
          stringsAsFactors = FALSE
        )
      )
    }

    data.frame(
      cohort = rep(
        cohort,
        length(missing_genes)
      ),
      gene = missing_genes,
      stringsAsFactors = FALSE
    )
  }
)

missing_signature <- do.call(
  rbind,
  missing_signature_list
)

if (is.null(missing_signature)) {
  missing_signature <- data.frame(
    cohort = character(0),
    gene = character(0),
    stringsAsFactors = FALSE
  )
}

write.csv(
  signature_coverage,
  file.path(
    RESULT_DIR,
    "09B_external_os_signature_coverage.csv"
  ),
  row.names = FALSE
)

write.csv(
  missing_signature,
  file.path(
    RESULT_DIR,
    "09B_external_os_missing_signature_genes.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Distribution audit after Terunuma gene collapse
# ------------------------------------------------------------

expression_distribution <- rbind(
  sample_distribution(
    tcga_ready$expression_complete,
    "TCGA"
  ),
  sample_distribution(
    teru_ready$expression_complete,
    "Terunuma_gene_median"
  ),
  sample_distribution(
    kao_ready$expression_complete,
    "Kao"
  )
)

write.csv(
  expression_distribution,
  file.path(
    RESULT_DIR,
    "09B_external_os_expression_distribution.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Save aligned complete-case objects for Step 10
# ------------------------------------------------------------

preprocessed <- list(
  settings = list(
    time_column = TIME_COL,
    event_column = EVENT_COL,
    age_column = AGE_COL,
    stage_column = STAGE_COL,
    Terunuma_mapping_rule = paste(
      "Remove probe IDs mapping to >1 unique symbol;",
      "collapse remaining probes per gene by sample-wise median."
    ),
    batch_correction_applied = FALSE
  ),
  TCGA = list(
    expression = tcga_ready$expression_complete,
    clinical = tcga_ready$clinical_complete
  ),
  Terunuma = list(
    expression = teru_ready$expression_complete,
    clinical = teru_ready$clinical_complete
  ),
  Kao = list(
    expression = kao_ready$expression_complete,
    clinical = kao_ready$clinical_complete
  ),
  common_genes = common_all_genes,
  expression_signature = expression_signature,
  mapping_audit = mapping_audit,
  sample_audit = sample_audit,
  gene_overlap = gene_overlap,
  signature_coverage = signature_coverage,
  expression_distribution = expression_distribution,
  created_at = Sys.time()
)

saveRDS(
  preprocessed,
  file.path(
    RESULT_DIR,
    "09B_external_os_preprocessed.rds"
  ),
  compress = "xz"
)


# ------------------------------------------------------------
# 9. Readiness
# ------------------------------------------------------------

common_signature_row <- signature_coverage[
  signature_coverage$cohort ==
    "Common_all",
  ,
  drop = FALSE
]

readiness <- data.frame(
  check = c(
    "Terunuma_all_input_IDs_present_in_annotation_DB",
    "Terunuma_unambiguous_gene_level_matrix_created",
    "TCGA_exact_overall_survival_columns_locked",
    "TCGA_has_complete_case_OS_cohort",
    "Terunuma_has_complete_case_OS_cohort",
    "Kao_has_complete_case_OS_cohort",
    "three_way_common_gene_set_nonempty",
    "119_gene_signature_loaded",
    "119_gene_signature_complete_in_all_three_cohorts",
    "batch_correction_applied"
  ),
  status = c(
    mapping_audit$value[
      mapping_audit$metric ==
        "IDs_present_in_annotation_DB"
    ] ==
      mapping_audit$value[
        mapping_audit$metric ==
          "input_transcript_clusters"
      ],
    nrow(teru_gene_exp) > 1000L,
    identical(
      TIME_COL,
      "OVERALL.SURVIVAL"
    ) &&
      identical(
        EVENT_COL,
        "overall.survival.indicator"
      ),
    ncol(
      tcga_ready$expression_complete
    ) > 0L,
    ncol(
      teru_ready$expression_complete
    ) > 0L,
    ncol(
      kao_ready$expression_complete
    ) > 0L,
    length(
      common_all_genes
    ) > 1000L,
    length(
      expression_signature
    ) == 119L,
    nrow(
      common_signature_row
    ) == 1L &&
      common_signature_row$genes_missing == 0L,
    FALSE
  ),
  note = c(
    "All transcript-cluster IDs exist as PROBEID keys.",
    paste0(
      nrow(teru_gene_exp),
      " gene symbols after removing ambiguous mappings and median collapse."
    ),
    paste(
      TIME_COL,
      "+",
      EVENT_COL
    ),
    paste0(
      ncol(
        tcga_ready$expression_complete
      ),
      " complete-case samples."
    ),
    paste0(
      ncol(
        teru_ready$expression_complete
      ),
      " complete-case samples."
    ),
    paste0(
      ncol(
        kao_ready$expression_complete
      ),
      " complete-case samples."
    ),
    paste0(
      length(
        common_all_genes
      ),
      " genes shared by all cohorts."
    ),
    paste0(
      length(
        expression_signature
      ),
      " genes loaded."
    ),
    paste0(
      common_signature_row$genes_missing,
      " signature genes missing from the three-way common set."
    ),
    "Intentionally FALSE. ComBat belongs to Step 10."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  readiness,
  file.path(
    RESULT_DIR,
    "09B_external_os_readiness.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Reproducibility report
# ------------------------------------------------------------

capture.output(
  {
    cat("STEP 09B — EXTERNAL OS PREPROCESSING\n")
    cat("Created: ", format(Sys.time()), "\n\n", sep = "")

    cat("MAPPING AUDIT\n")
    print(mapping_audit)

    cat("\nSAMPLE AUDIT\n")
    print(sample_audit)

    cat("\nGENE OVERLAP\n")
    print(gene_overlap)

    cat("\nSIGNATURE COVERAGE\n")
    print(signature_coverage)

    cat("\nEXPRESSION DISTRIBUTION\n")
    print(expression_distribution)

    cat("\nREADINESS\n")
    print(readiness)

    cat("\nSESSION INFO\n")
    sessionInfo()
  },
  file = file.path(
    RESULT_DIR,
    "09B_external_os_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 11. Console output
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 09B COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nTERUNUMA MAPPING AUDIT\n")
print(
  mapping_audit,
  row.names = FALSE
)

cat("\nOS SAMPLE AUDIT — EXACT OVERALL.SURVIVAL\n")
print(
  sample_audit,
  row.names = FALSE
)

cat("\nGENE OVERLAP AFTER TERUNUMA MEDIAN COLLAPSE\n")
print(
  gene_overlap,
  row.names = FALSE
)

cat("\n119-GENE SIGNATURE COVERAGE\n")
print(
  signature_coverage,
  row.names = FALSE
)

cat("\nEXPRESSION DISTRIBUTIONS\n")
print(
  expression_distribution,
  row.names = FALSE
)

cat("\nREADINESS\n")
print(
  readiness,
  row.names = FALSE
)

cat(
  "\nPreprocessed object for Step 10:\n",
  file.path(
    RESULT_DIR,
    "09B_external_os_preprocessed.rds"
  ),
  "\n",
  sep = ""
)

cat("\nFinished successfully.\n")
