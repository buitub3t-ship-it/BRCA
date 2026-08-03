#!/usr/bin/env Rscript

# ============================================================
# RUN ORIGINAL M2EFM get_m2eqtls() CODE
# Run from repository root:
# Rscript run_get_m2eqtls_author.R
# ============================================================


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

ROOT_DIR <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

R_DIR <- file.path(ROOT_DIR, "R")

DATA_DIR <- file.path(
  ROOT_DIR,
  "inst",
  "extdata",
  "csv_output"
)

RESULT_DIR <- file.path(ROOT_DIR, "results")

dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Repository root:", ROOT_DIR, "\n")


# ------------------------------------------------------------
# 2. CHECK INPUT FILES
# ------------------------------------------------------------

EXP_FILE <- file.path(
  DATA_DIR,
  "TCGA_BRCA_EXP.csv"
)

METH_FILE <- file.path(
  DATA_DIR,
  "TCGA_BRCA_METH.csv"
)

PROBE_CSV <- file.path(
  DATA_DIR,
  "PROBELIST.csv"
)

required_files <- c(
  EXP_FILE,
  METH_FILE,
  PROBE_CSV,
  file.path(R_DIR, "package_loader.R"),
  file.path(R_DIR, "ProfileData.R"),
  file.path(R_DIR, "MethylationData.R"),
  file.path(R_DIR, "ExpressionData.R"),
  file.path(R_DIR, "eqtl.R"),
  file.path(R_DIR, "m2eQTL.R")
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing files:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------
# 3. CHECK PACKAGES
#
# All packages must already be installed so that the author's
# old package_loader.R does not call biocLite().
# ------------------------------------------------------------

required_packages <- c(
  "R6",
  "data.table",
  "impute",
  "lumi",
  "MatrixEQTL",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "minfi",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

installed <- vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

if (!all(installed)) {
  stop(
    "Missing R packages:\n",
    paste(names(installed)[!installed], collapse = "\n"),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(impute)
  library(lumi)
  library(MatrixEQTL)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})


# ------------------------------------------------------------
# 4. SOURCE ORIGINAL AUTHOR CODE
#
# Không chỉnh sửa lại thuật toán trong các file này.
# ------------------------------------------------------------

source(
  file.path(R_DIR, "package_loader.R"),
  local = .GlobalEnv
)

source(
  file.path(R_DIR, "ProfileData.R"),
  local = .GlobalEnv
)

source(
  file.path(R_DIR, "MethylationData.R"),
  local = .GlobalEnv
)

source(
  file.path(R_DIR, "ExpressionData.R"),
  local = .GlobalEnv
)

source(
  file.path(R_DIR, "eqtl.R"),
  local = .GlobalEnv
)

source(
  file.path(R_DIR, "m2eQTL.R"),
  local = .GlobalEnv
)

required_objects <- c(
  "ProfileData",
  "MethylationData",
  "ExpressionData",
  "get_eqtls",
  "get_m2eqtls"
)

object_status <- vapply(
  required_objects,
  exists,
  logical(1),
  inherits = TRUE
)

if (!all(object_status)) {
  stop(
    "Objects not loaded:\n",
    paste(names(object_status)[!object_status], collapse = "\n")
  )
}

cat("Original author functions loaded successfully.\n")


# ------------------------------------------------------------
# 5. CONVERT PROBELIST TO AUTHOR FORMAT
#
# get_m2eqtls() uses:
# read.table(probe_list_loc, header = FALSE)
#
# Therefore create one-column TXT without a header.
# Order is preserved.
# ------------------------------------------------------------

probe_table <- data.table::fread(
  PROBE_CSV,
  data.table = FALSE,
  check.names = FALSE
)

probe_ids <- trimws(
  as.character(probe_table[[1L]])
)

probe_ids <- probe_ids[
  !is.na(probe_ids) &
  nzchar(probe_ids)
]

if (anyDuplicated(probe_ids)) {
  stop("Duplicated probes found in PROBELIST.csv.")
}

PROBE_TXT <- file.path(
  RESULT_DIR,
  "PROBELIST_author_format.txt"
)

write.table(
  probe_ids,
  file = PROBE_TXT,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("Probe-list entries:", length(probe_ids), "\n")


# ------------------------------------------------------------
# 6. RUN ORIGINAL get_m2eqtls()
#
# Paper/vignette settings:
#   550 candidate probes
#   top 110 trans associations
#
# Other settings match the original function.
# ------------------------------------------------------------

start_time <- Sys.time()

cat("Starting get_m2eqtls at:", format(start_time), "\n")

m2e_os <- get_m2eqtls(
  probe_list_loc = PROBE_TXT,
  meth_data       = METH_FILE,
  exp_data        = EXP_FILE,
  gene_list       = NULL,

  num_cis         = 1e10,
  num_trans       = 110,
  num_probes      = 550,

  treatment_col   = "-01",
  cc              = 1e-3,
  tc              = 1e-10
)

end_time <- Sys.time()

elapsed_minutes <- as.numeric(
  difftime(
    end_time,
    start_time,
    units = "mins"
  )
)

cat(
  "get_m2eqtls completed in",
  round(elapsed_minutes, 2),
  "minutes.\n"
)


# ------------------------------------------------------------
# 7. BASIC OUTPUT CHECK
# ------------------------------------------------------------

if (
  !is.list(m2e_os) ||
  !all(c("result", "cis_link", "trans_link") %in% names(m2e_os))
) {
  stop("Unexpected get_m2eqtls output structure.")
}

n_cis <- nrow(m2e_os$result$cis$eqtls)
n_trans <- nrow(m2e_os$result$trans$eqtls)
n_cis_links <- nrow(m2e_os$cis_link)
n_trans_links <- nrow(m2e_os$trans_link)

summary_table <- data.frame(
  metric = c(
    "cis_eqtls",
    "trans_eqtls",
    "cis_links",
    "trans_links",
    "elapsed_minutes"
  ),
  value = c(
    n_cis,
    n_trans,
    n_cis_links,
    n_trans_links,
    elapsed_minutes
  )
)

print(summary_table, row.names = FALSE)


# ------------------------------------------------------------
# 8. SAVE RESULTS
# ------------------------------------------------------------

saveRDS(
  m2e_os,
  file.path(
    RESULT_DIR,
    "01_TCGA_OS_m2eQTL.rds"
  )
)

write.csv(
  m2e_os$result$cis$eqtls,
  file.path(
    RESULT_DIR,
    "01_cis_m2eQTL.csv"
  ),
  row.names = FALSE
)

write.csv(
  m2e_os$result$trans$eqtls,
  file.path(
    RESULT_DIR,
    "01_trans_m2eQTL.csv"
  ),
  row.names = FALSE
)

write.csv(
  m2e_os$cis_link,
  file.path(
    RESULT_DIR,
    "01_cis_links.csv"
  ),
  row.names = FALSE
)

write.csv(
  m2e_os$trans_link,
  file.path(
    RESULT_DIR,
    "01_trans_links.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_table,
  file.path(
    RESULT_DIR,
    "01_m2eQTL_summary.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    RESULT_DIR,
    "01_sessionInfo.txt"
  )
)

cat(
  "Results saved to:",
  normalizePath(RESULT_DIR, winslash = "/"),
  "\n"
)
