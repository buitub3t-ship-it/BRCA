# Run from the M2EFM_OS_TCGA_KAO_TERU project root.
# Package installation is intentionally separate: run 00_install_packages.R once.

scripts <- c(
  "01_preflight_and_prepare.R",
  "02_identify_m2eqtls.R",
  "03_train_internal_full.R",
  "04_train_external_expression.R",
  "05_save_session_info.R"
)

for (script in scripts) {
  message("\n==============================")
  message("Running ", script)
  message("==============================")
  source(script, local = new.env(parent = globalenv()))
  gc()
}

message("All OS-M2EFM steps completed. See ./results.")
