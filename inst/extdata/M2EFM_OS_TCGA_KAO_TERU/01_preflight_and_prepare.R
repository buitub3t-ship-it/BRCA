source("config.R")
source(file.path("R", "io_and_preprocessing.R"))

assert_packages(c("data.table", "AnnotationDbi", "hugene10sttranscriptcluster.db"))
message("Data directory: ", DATA_DIR)

# Load raw matrices and clinical data.
tcga_exp <- read_feature_matrix(FILES$tcga_exp, "gene_id")
tcga_meth <- read_feature_matrix(FILES$tcga_meth, "probe_id")
kao_exp <- read_feature_matrix(FILES$kao_exp, "gene_id")
teru_probe_exp <- read_feature_matrix(FILES$teru_exp, "probe_id")

tcga_clin <- clean_clinical(read_clinical(FILES$tcga_clin), "TCGA")
kao_clin <- clean_clinical(read_clinical(FILES$kao_clin), "KAO")
teru_clin <- clean_clinical(read_clinical(FILES$teru_clin), "TERU")

probelist_df <- data.table::fread(FILES$probelist, data.table = FALSE)
if (!"probe_id" %in% names(probelist_df)) {
  stop("PROBELIST.csv must contain column probe_id.")
}
probelist <- unique(as.character(probelist_df$probe_id))

if (length(probelist) < NUM_DISCOVERY_PROBES) {
  stop("PROBELIST has fewer than ", NUM_DISCOVERY_PROBES, " probes.")
}
missing_top <- setdiff(head(probelist, NUM_DISCOVERY_PROBES), rownames(tcga_meth))
if (length(missing_top) > 0L) {
  stop("Top 550 probes missing from TCGA methylation: ",
       paste(missing_top, collapse = ", "))
}

message("Mapping Terunuma Human Gene 1.0 ST transcript clusters to gene symbols...")
teru_mapping <- map_terunuma_to_gene(teru_probe_exp)
teru_exp <- teru_mapping$expression

# The paper filtered TCGA genes to those available in the validation data.
common_genes <- Reduce(
  intersect,
  list(rownames(tcga_exp), rownames(kao_exp), rownames(teru_exp))
)
if (length(common_genes) < 1000L) {
  stop("Only ", length(common_genes), " common genes remain across cohorts.")
}

# Preserve TCGA row order.
common_genes <- rownames(tcga_exp)[rownames(tcga_exp) %in% common_genes]
tcga_exp <- tcga_exp[common_genes, , drop = FALSE]
kao_exp <- kao_exp[common_genes, , drop = FALSE]
teru_exp <- teru_exp[common_genes, , drop = FALSE]

# Align molecular and clinical sample IDs.
tcga_aligned <- align_clinical_to_expression(tcga_clin, tcga_exp)
tcga_exp <- tcga_aligned$expression
tcga_clin <- tcga_aligned$clinical

kao_aligned <- align_clinical_to_expression(kao_clin, kao_exp)
kao_exp <- kao_aligned$expression
kao_clin <- kao_aligned$clinical

teru_aligned <- align_clinical_to_expression(teru_clin, teru_exp)
teru_exp <- teru_aligned$expression
teru_clin <- teru_aligned$clinical

if (APPLY_COMBAT_TO_CURRENT_MATRICES) {
  stop(
    "APPLY_COMBAT_TO_CURRENT_MATRICES=TRUE is intentionally blocked. ",
    "The supplied TCGA/Kao matrices appear already harmonized; re-ComBat ",
    "would double-correct them. Rebuild all three cohorts from pre-ComBat ",
    "expression matrices before enabling this step."
  )
}

# Current uploaded-file audit.
audit <- rbind(
  make_audit_row("TCGA_expression", tcga_exp, tcga_clin),
  make_audit_row("TCGA_methylation", tcga_meth, tcga_clin),
  make_audit_row("Kao_expression", kao_exp, kao_clin),
  make_audit_row("Terunuma_expression_gene_level", teru_exp, teru_clin)
)
audit$notes <- c(
  paste0("Common-gene matrix; complete OS expected near 1039 in current files."),
  paste0("Contains tumor/normal/metastatic samples; matched tumor subset selected later."),
  "External validation cohort.",
  paste0(
    "Mapped from ", teru_mapping$n_input_probes, " transcript clusters; ",
    teru_mapping$n_gene_symbols, " gene symbols before common-gene filtering."
  )
)

if (STRICT_PAPER_COUNTS) {
  n_tcga_complete <- sum(complete_os_rows(tcga_clin))
  if (nrow(tcga_exp) != 10990L || n_tcga_complete != 1028L) {
    stop(
      "STRICT_PAPER_COUNTS is TRUE, but current data have ", nrow(tcga_exp),
      " common genes and ", n_tcga_complete,
      " TCGA expression samples with complete OS."
    )
  }
}

saveRDS(tcga_exp, file.path(CACHE_DIR, "tcga_exp_common.rds"), compress = FALSE)
saveRDS(tcga_meth, file.path(CACHE_DIR, "tcga_meth_beta.rds"), compress = FALSE)
saveRDS(kao_exp, file.path(CACHE_DIR, "kao_exp_common.rds"), compress = FALSE)
saveRDS(teru_exp, file.path(CACHE_DIR, "teru_exp_common.rds"), compress = FALSE)
saveRDS(tcga_clin, file.path(CACHE_DIR, "tcga_clin_clean.rds"))
saveRDS(kao_clin, file.path(CACHE_DIR, "kao_clin_clean.rds"))
saveRDS(teru_clin, file.path(CACHE_DIR, "teru_clin_clean.rds"))
saveRDS(probelist, file.path(CACHE_DIR, "probelist.rds"))
saveRDS(teru_mapping$mapping, file.path(CACHE_DIR, "teru_probe_to_symbol.rds"))

data.table::fwrite(audit, file.path(RESULTS_DIR, "01_data_audit.csv"))
data.table::fwrite(
  data.frame(gene_id = common_genes),
  file.path(RESULTS_DIR, "01_common_genes.csv")
)

data.table::fwrite(
  data.frame(
    item = c(
      "PROBELIST_total",
      "PROBELIST_present_in_methylation",
      "PROBELIST_missing_total",
      "top_550_missing",
      "common_genes_TCGA_Kao_Terunuma"
    ),
    value = c(
      length(probelist),
      sum(probelist %in% rownames(tcga_meth)),
      sum(!probelist %in% rownames(tcga_meth)),
      length(missing_top),
      length(common_genes)
    )
  ),
  file.path(RESULTS_DIR, "01_preflight_summary.csv")
)

message("Preflight/preparation completed.")
print(audit)
