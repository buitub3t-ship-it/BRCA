source("config.R")
source(file.path("R", "io_and_preprocessing.R"))
source(file.path("R", "m2eqtl.R"))

required_cache <- c(
  "tcga_exp_common.rds", "tcga_meth_beta.rds", "probelist.rds"
)
if (!all(file.exists(file.path(CACHE_DIR, required_cache)))) {
  stop("Run 01_preflight_and_prepare.R first.")
}

tcga_exp <- readRDS(file.path(CACHE_DIR, "tcga_exp_common.rds"))
tcga_meth <- readRDS(file.path(CACHE_DIR, "tcga_meth_beta.rds"))
probelist <- readRDS(file.path(CACHE_DIR, "probelist.rds"))

m2e <- run_m2eqtl_discovery(
  tcga_expression = tcga_exp,
  tcga_methylation_beta = tcga_meth,
  probelist = probelist
)

saveRDS(m2e, file.path(CACHE_DIR, "m2eqtl_discovery.rds"), compress = TRUE)

data.table::fwrite(
  m2e$trans_selected,
  file.path(RESULTS_DIR, "02_trans_m2eqtl_selected.csv")
)
data.table::fwrite(
  m2e$cis_replacements,
  file.path(RESULTS_DIR, "02_cis_replacement_m2eqtl.csv")
)
data.table::fwrite(
  data.frame(gene_id = m2e$trans_genes),
  file.path(RESULTS_DIR, "02_trans_m2egenes.csv")
)
data.table::fwrite(
  data.frame(gene_id = m2e$cis_genes),
  file.path(RESULTS_DIR, "02_cis_m2egenes.csv")
)
data.table::fwrite(
  data.frame(gene_id = m2e$expression_signature),
  file.path(RESULTS_DIR, "02_expression_only_signature.csv")
)
data.table::fwrite(
  data.frame(probe_id = m2e$methylation_probes),
  file.path(RESULTS_DIR, "02_m2eqtl_probes.csv")
)
data.table::fwrite(
  data.frame(sample_id = m2e$matched_samples),
  file.path(RESULTS_DIR, "02_matched_TCGA_discovery_samples.csv")
)

summary_table <- data.frame(
  item = c(
    "matched_TCGA_tumor_samples",
    "discovery_CpGs",
    "significant_trans_pairs",
    "selected_trans_genes",
    "selected_m2eQTL_probes",
    "cis_replacement_genes",
    "expression_only_signature_size"
  ),
  value = c(
    length(m2e$matched_samples),
    length(m2e$discovery_probes),
    nrow(m2e$trans_all_significant),
    length(m2e$trans_genes),
    length(m2e$methylation_probes),
    length(m2e$cis_genes),
    length(m2e$expression_signature)
  )
)
data.table::fwrite(
  summary_table,
  file.path(RESULTS_DIR, "02_m2eqtl_summary.csv")
)

message("m2eQTL discovery completed.")
print(summary_table)
