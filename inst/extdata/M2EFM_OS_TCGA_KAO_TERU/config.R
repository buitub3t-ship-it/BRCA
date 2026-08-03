# Configuration for OS-M2EFM replication using TCGA, Kao, and Terunuma.
# Put this project folder next to csv_output/, or set M2EFM_DATA_DIR.

options(stringsAsFactors = FALSE)

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

required_files <- c(
  "TCGA_BRCA_EXP.csv",
  "TCGA_BRCA_METH.csv",
  "TCGA_BRCA_CLIN.csv",
  "KAO_BRCA_EXP.csv",
  "KAO_BRCA_CLIN.csv",
  "TERUNUMA_BRCA_EXP.csv",
  "TERUNUMA_BRCA_CLIN.csv",
  "PROBELIST.csv"
)

resolve_data_dir <- function() {
  env_dir <- Sys.getenv("M2EFM_DATA_DIR", unset = "")
  candidates <- unique(c(
    env_dir,
    file.path(PROJECT_ROOT, "csv_output"),
    PROJECT_ROOT
  ))
  candidates <- candidates[nzchar(candidates)]

  for (candidate in candidates) {
    if (dir.exists(candidate) &&
        all(file.exists(file.path(candidate, required_files)))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  stop(
    "Cannot locate all input CSV files. Put them in ./csv_output, ",
    "in the project root, or set environment variable M2EFM_DATA_DIR."
  )
}

DATA_DIR <- resolve_data_dir()
CACHE_DIR <- file.path(PROJECT_ROOT, "cache")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

FILES <- list(
  tcga_exp = file.path(DATA_DIR, "TCGA_BRCA_EXP.csv"),
  tcga_meth = file.path(DATA_DIR, "TCGA_BRCA_METH.csv"),
  tcga_clin = file.path(DATA_DIR, "TCGA_BRCA_CLIN.csv"),
  kao_exp = file.path(DATA_DIR, "KAO_BRCA_EXP.csv"),
  kao_clin = file.path(DATA_DIR, "KAO_BRCA_CLIN.csv"),
  teru_exp = file.path(DATA_DIR, "TERUNUMA_BRCA_EXP.csv"),
  teru_clin = file.path(DATA_DIR, "TERUNUMA_BRCA_CLIN.csv"),
  probelist = file.path(DATA_DIR, "PROBELIST.csv")
)

# Paper-oriented settings.
NUM_DISCOVERY_PROBES <- 550L
NUM_TRANS_GENES <- 110L
CIS_P_THRESHOLD <- 1e-3
TRANS_P_THRESHOLD <- 1e-10
CIS_DISTANCE_BP <- 10000L
MAD_THRESHOLD <- 0.05

# Model evaluation settings.
N_MONTE_CARLO_SPLITS <- 100L
TRAINING_PROPORTION <- 0.70
GLMNET_CV_FOLDS <- 10L
BASE_SEED <- 1L

# TRUE reproduces the author's evaluate() behavior: calculate per-feature
# mean/range once on all TCGA samples before the 70/30 Monte-Carlo splits.
AUTHOR_STYLE_GLOBAL_SCALING <- TRUE

# Current TCGA and Kao matrices appear to be the author's already harmonized
# matrices. Re-running ComBat would double-correct them, so leave this FALSE.
APPLY_COMBAT_TO_CURRENT_MATRICES <- FALSE

# The uploaded repo-style data contain 10,000 genes and 1,039 TCGA expression
# samples with complete OS, not the paper's 10,990 genes / 1,028 samples.
# Set TRUE only when you have an exact paper cohort/sample list and matrices.
STRICT_PAPER_COUNTS <- FALSE

STAGE_LEVELS <- c("Stage I", "Stage II", "Stage III", "Stage IV")
