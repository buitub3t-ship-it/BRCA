#!/usr/bin/env Rscript

# ============================================================
# STEP 10B — PAM50 rorS FOR FIGURE 2B
#
# Input:
#   results/10_external_os_combat_ready.rds
#
# Primary implementation:
#   Reproduce the repository helper get_pam50() up to the rorS call:
#     - map gene aliases to Entrez IDs with org.Hs.egALIAS2EG;
#     - transpose expression to sample x gene;
#     - call genefu::rorS(data = ..., annot = ...).
#
# rorS is computed separately for TCGA, Terunuma, and Kao, matching
# the repository helper's dataset-at-a-time design. A combined-cohort
# calculation is retained only as a sensitivity audit.
#
# This script does not call molecular.subtyping(), because that optional
# function is not required for the rorS comparison and previously failed
# in this installed genefu environment due to unavailable pam50.robust.
#
# Outputs:
#   results/10B_external_os_rorS.csv
#   results/10B_external_os_rorS_audit.csv
#   results/10B_external_os_rorS_mapping_audit.csv
#   results/10B_external_os_rorS_result.rds
#   results/10B_external_os_combat_rors_ready.rds
#   results/10B_external_os_readiness.csv
#   results/10B_external_os_sessionInfo.txt
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

input_file <- file.path(
  result_dir,
  "10_external_os_combat_ready.rds"
)

required_packages <- c(
  "genefu",
  "AnnotationDbi",
  "org.Hs.eg.db"
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
  library(genefu)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

if (!file.exists(input_file)) {
  stop(
    "Missing input: ",
    input_file
  )
}


# ------------------------------------------------------------
# 1. Helpers
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
      " expression contains duplicated gene/sample names."
    )
  }

  if (
    anyNA(x) ||
      any(!is.finite(x))
  ) {
    stop(
      cohort,
      " expression contains missing/non-finite values."
    )
  }

  x
}

build_author_style_annotation <- function(
    genes
) {

  alias_map <- org.Hs.eg.db::org.Hs.egALIAS2EG

  mapped_aliases <- AnnotationDbi::mappedkeys(
    alias_map
  )

  alias_list <- AnnotationDbi::as.list(
    alias_map[
      mapped_aliases
    ]
  )

  mapped_genes <- genes[
    genes %in%
      names(alias_list)
  ]

  entrez_ids <- vapply(
    alias_list[
      mapped_genes
    ],
    function(values) {
      as.character(
        values[[1L]]
      )
    },
    character(1)
  )

  annotation <- data.frame(
    GeneSymbol = mapped_genes,
    EntrezGene.ID = unname(
      entrez_ids
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(annotation) <- mapped_genes

  annotation
}

extract_score <- function(
    result,
    sample_ids
) {

  if (
    !is.list(result) ||
      is.null(result$score)
  ) {
    stop(
      "genefu::rorS did not return $score."
    )
  }

  score_object <- result$score

  if (
    is.data.frame(score_object) ||
      is.matrix(score_object)
  ) {

    if (ncol(score_object) < 1L) {
      stop(
        "rorS score object has zero columns."
      )
    }

    score_names <- rownames(
      score_object
    )

    score <- as.numeric(
      score_object[
        ,
        1L
      ]
    )

    if (
      !is.null(score_names) &&
        length(score_names) ==
          length(score)
    ) {
      names(score) <- score_names
    }

  } else {

    score <- as.numeric(
      score_object
    )

    if (
      !is.null(
        names(score_object)
      )
    ) {
      names(score) <- names(
        score_object
      )
    }
  }

  if (length(score) != length(sample_ids)) {
    stop(
      "rorS returned ",
      length(score),
      " scores for ",
      length(sample_ids),
      " samples."
    )
  }

  if (
    is.null(names(score)) ||
      anyNA(names(score)) ||
      !setequal(
        names(score),
        sample_ids
      )
  ) {
    names(score) <- sample_ids
  } else {
    score <- score[
      sample_ids
    ]
  }

  if (
    anyNA(score) ||
      any(!is.finite(score))
  ) {
    stop(
      "rorS scores contain missing/non-finite values."
    )
  }

  score
}

compute_rors <- function(
    expression,
    annotation,
    label
) {

  shared_genes <- intersect(
    rownames(expression),
    rownames(annotation)
  )

  expression_mapped <- expression[
    shared_genes,
    ,
    drop = FALSE
  ]

  annotation_mapped <- annotation[
    shared_genes,
    ,
    drop = FALSE
  ]

  if (!identical(
    rownames(expression_mapped),
    rownames(annotation_mapped)
  )) {
    stop(
      "Expression/annotation gene order mismatch for ",
      label,
      "."
    )
  }

  sample_by_gene <- t(
    expression_mapped
  )

  cat(
    "Computing rorS for ",
    label,
    ": ",
    nrow(sample_by_gene),
    " samples x ",
    ncol(sample_by_gene),
    " mapped genes...\n",
    sep = ""
  )

  result <- tryCatch(
    genefu::rorS(
      data = sample_by_gene,
      annot = annotation_mapped
    ),
    error = function(error_condition) {
      stop(
        "rorS failed for ",
        label,
        ": ",
        conditionMessage(
          error_condition
        )
      )
    }
  )

  score <- extract_score(
    result,
    rownames(sample_by_gene)
  )

  list(
    result = result,
    score = score,
    genes_used = shared_genes
  )
}


# ------------------------------------------------------------
# 2. Load and validate Step 10A object
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

cohorts <- c(
  "TCGA",
  "Terunuma",
  "Kao"
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

    x <- ready[[cohort]]$clinical

    if (is.null(rownames(x))) {
      stop(
        cohort,
        " clinical data lacks sample row names."
      )
    }

    if (!identical(
      colnames(expression[[cohort]]),
      rownames(x)
    )) {
      stop(
        cohort,
        " expression and clinical samples are not aligned."
      )
    }

    x
  }
)

names(clinical) <- cohorts

all_genes <- rownames(
  expression$TCGA
)

for (cohort in cohorts) {
  if (!identical(
    rownames(expression[[cohort]]),
    all_genes
  )) {
    stop(
      "Gene order differs for ",
      cohort,
      "."
    )
  }
}


# ------------------------------------------------------------
# 3. Build the repository-style alias-to-Entrez annotation
# ------------------------------------------------------------

annotation <- build_author_style_annotation(
  all_genes
)

mapping_audit <- data.frame(
  metric = c(
    "ComBat_common_genes",
    "genes_mapped_by_org.Hs.egALIAS2EG",
    "genes_not_mapped",
    "mapping_fraction",
    "mapping_database",
    "mapping_database_version",
    "genefu_version"
  ),
  value = c(
    length(all_genes),
    nrow(annotation),
    length(
      setdiff(
        all_genes,
        rownames(annotation)
      )
    ),
    nrow(annotation) /
      length(all_genes),
    "org.Hs.eg.db::org.Hs.egALIAS2EG",
    as.character(
      packageVersion(
        "org.Hs.eg.db"
      )
    ),
    as.character(
      packageVersion(
        "genefu"
      )
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  annotation,
  file.path(
    result_dir,
    "10B_external_os_rorS_gene_annotation.csv"
  ),
  row.names = FALSE
)

write.csv(
  mapping_audit,
  file.path(
    result_dir,
    "10B_external_os_rorS_mapping_audit.csv"
  ),
  row.names = FALSE
)

writeLines(
  setdiff(
    all_genes,
    rownames(annotation)
  ),
  file.path(
    result_dir,
    "10B_external_os_rorS_unmapped_genes.txt"
  )
)

if (nrow(annotation) < 100L) {
  stop(
    "Too few genes mapped for rorS: ",
    nrow(annotation),
    "."
  )
}


# ------------------------------------------------------------
# 4. Primary author-style rorS: one call per dataset
# ------------------------------------------------------------

primary_results <- lapply(
  cohorts,
  function(cohort) {
    compute_rors(
      expression[[cohort]],
      annotation,
      cohort
    )
  }
)

names(primary_results) <- cohorts


# ------------------------------------------------------------
# 5. Sensitivity rorS: one call on the combined corrected matrix
# ------------------------------------------------------------

combined_expression <- do.call(
  cbind,
  expression
)

if (anyDuplicated(
  colnames(combined_expression)
)) {
  stop(
    "Duplicated sample IDs across cohorts."
  )
}

combined_result <- compute_rors(
  combined_expression,
  annotation,
  "combined TCGA + Terunuma + Kao"
)


# ------------------------------------------------------------
# 6. Create score table and attach primary rorS to clinical data
# ------------------------------------------------------------

score_rows <- list()

for (cohort in cohorts) {

  sample_ids <- colnames(
    expression[[cohort]]
  )

  primary_score <- primary_results[[cohort]]$score[
    sample_ids
  ]

  combined_score <- combined_result$score[
    sample_ids
  ]

  if (
    anyNA(primary_score) ||
      anyNA(combined_score)
  ) {
    stop(
      "Failed to align rorS scores for ",
      cohort,
      "."
    )
  }

  score_rows[[cohort]] <- data.frame(
    sample_id = sample_ids,
    cohort = cohort,
    rorS = as.numeric(
      primary_score
    ),
    rorS_combined_sensitivity = as.numeric(
      combined_score
    ),
    stringsAsFactors = FALSE
  )

  clinical[[cohort]]$rorS <- as.numeric(
    primary_score[
      rownames(
        clinical[[cohort]]
      )
    ]
  )

  if (
    anyNA(
      clinical[[cohort]]$rorS
    )
  ) {
    stop(
      "Failed to attach rorS to clinical data for ",
      cohort,
      "."
    )
  }
}

score_table <- do.call(
  rbind,
  score_rows
)

rownames(score_table) <- NULL

write.csv(
  score_table,
  file.path(
    result_dir,
    "10B_external_os_rorS.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Audit distributions and sensitivity agreement
# ------------------------------------------------------------

audit <- do.call(
  rbind,
  lapply(
    cohorts,
    function(cohort) {

      x <- score_table[
        score_table$cohort == cohort,
        ,
        drop = FALSE
      ]

      data.frame(
        cohort = cohort,
        samples = nrow(x),
        rorS_min = min(x$rorS),
        rorS_q25 = unname(
          quantile(
            x$rorS,
            0.25
          )
        ),
        rorS_median = median(
          x$rorS
        ),
        rorS_mean = mean(
          x$rorS
        ),
        rorS_q75 = unname(
          quantile(
            x$rorS,
            0.75
          )
        ),
        rorS_max = max(
          x$rorS
        ),
        rorS_sd = stats::sd(
          x$rorS
        ),
        combined_sensitivity_min = min(
          x$rorS_combined_sensitivity
        ),
        combined_sensitivity_median = median(
          x$rorS_combined_sensitivity
        ),
        combined_sensitivity_mean = mean(
          x$rorS_combined_sensitivity
        ),
        combined_sensitivity_max = max(
          x$rorS_combined_sensitivity
        ),
        primary_vs_combined_spearman = suppressWarnings(
          stats::cor(
            x$rorS,
            x$rorS_combined_sensitivity,
            method = "spearman"
          )
        ),
        primary_vs_combined_pearson = suppressWarnings(
          stats::cor(
            x$rorS,
            x$rorS_combined_sensitivity,
            method = "pearson"
          )
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  audit,
  file.path(
    result_dir,
    "10B_external_os_rorS_audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Save result objects
# ------------------------------------------------------------

result_object <- list(
  primary_per_cohort = lapply(
    primary_results,
    function(x) {
      x$result
    }
  ),
  combined_sensitivity = combined_result$result,
  scores = score_table,
  annotation = annotation,
  mapping_audit = mapping_audit,
  audit = audit,
  implementation = paste(
    "Primary rorS computed separately per cohort using",
    "org.Hs.egALIAS2EG mapping and genefu::rorS;",
    "combined-cohort result retained only as sensitivity."
  ),
  created_at = Sys.time()
)

saveRDS(
  result_object,
  file.path(
    result_dir,
    "10B_external_os_rorS_result.rds"
  ),
  compress = "xz"
)

ready$TCGA$clinical <- clinical$TCGA
ready$Terunuma$clinical <- clinical$Terunuma
ready$Kao$clinical <- clinical$Kao

ready$rorS <- score_table
ready$rorS_mapping_audit <- mapping_audit
ready$rorS_audit <- audit
ready$settings$rorS_computed <- TRUE
ready$settings$rorS_primary_method <- paste(
  "Separate genefu::rorS call per cohort after joint ComBat,",
  "using org.Hs.egALIAS2EG, matching repository get_pam50() logic."
)
ready$settings$rorS_combined_sensitivity_available <- TRUE
ready$created_at_10B <- Sys.time()

output_ready <- file.path(
  result_dir,
  "10B_external_os_combat_rors_ready.rds"
)

saveRDS(
  ready,
  output_ready,
  compress = "xz"
)


# ------------------------------------------------------------
# 9. Readiness
# ------------------------------------------------------------

readiness <- data.frame(
  check = c(
    "Step_10A_object_loaded",
    "common_gene_order_preserved",
    "author_style_alias_mapping_created",
    "rorS_TCGA_complete",
    "rorS_Terunuma_complete",
    "rorS_Kao_complete",
    "rorS_combined_sensitivity_complete",
    "rorS_attached_to_all_clinical_tables",
    "Step_11_ready_object_saved",
    "NCA_signature_available"
  ),
  status = c(
    TRUE,
    TRUE,
    nrow(annotation) >= 100L,
    sum(
      score_table$cohort == "TCGA"
    ) == ncol(
      expression$TCGA
    ) &&
      all(
        is.finite(
          score_table$rorS[
            score_table$cohort == "TCGA"
          ]
        )
      ),
    sum(
      score_table$cohort == "Terunuma"
    ) == ncol(
      expression$Terunuma
    ) &&
      all(
        is.finite(
          score_table$rorS[
            score_table$cohort == "Terunuma"
          ]
        )
      ),
    sum(
      score_table$cohort == "Kao"
    ) == ncol(
      expression$Kao
    ) &&
      all(
        is.finite(
          score_table$rorS[
            score_table$cohort == "Kao"
          ]
        )
      ),
    all(
      is.finite(
        score_table$rorS_combined_sensitivity
      )
    ),
    all(
      vapply(
        clinical,
        function(x) {
          "rorS" %in%
            colnames(x) &&
            all(
              is.finite(
                x$rorS
              )
            )
        },
        logical(1)
      )
    ),
    file.exists(output_ready),
    FALSE
  ),
  note = c(
    normalizePath(
      input_file,
      winslash = "/",
      mustWork = TRUE
    ),
    paste0(
      length(all_genes),
      " genes."
    ),
    paste0(
      nrow(annotation),
      "/",
      length(all_genes),
      " genes mapped."
    ),
    paste0(
      ncol(expression$TCGA),
      " scores."
    ),
    paste0(
      ncol(expression$Terunuma),
      " scores."
    ),
    paste0(
      ncol(expression$Kao),
      " scores."
    ),
    paste0(
      nrow(score_table),
      " sensitivity scores."
    ),
    "Primary per-cohort rorS aligned by sample ID.",
    output_ready,
    "FALSE: exact NCA signature remains unavailable."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  readiness,
  file.path(
    result_dir,
    "10B_external_os_readiness.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Session report and console summary
# ------------------------------------------------------------

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
    cat("STEP 10B — PAM50 rorS\n")
    cat("Created: ", format(end_time), "\n", sep = "")
    cat("Elapsed minutes: ", elapsed_minutes, "\n\n", sep = "")

    cat("MAPPING AUDIT\n")
    print(mapping_audit)

    cat("\nrorS AUDIT\n")
    print(audit)

    cat("\nREADINESS\n")
    print(readiness)

    cat("\nSESSION INFO\n")
    sessionInfo()
  },
  file = file.path(
    result_dir,
    "10B_external_os_sessionInfo.txt"
  )
)

cat("\n")
cat("============================================================\n")
cat("STEP 10B COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nMAPPING AUDIT\n")
print(
  mapping_audit,
  row.names = FALSE
)

cat("\nrorS AUDIT\n")
print(
  audit,
  row.names = FALSE
)

cat("\nREADINESS\n")
print(
  readiness,
  row.names = FALSE
)

cat(
  "\nStep 11 ready object:\n",
  output_ready,
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

cat("\nFinished successfully.\n")
