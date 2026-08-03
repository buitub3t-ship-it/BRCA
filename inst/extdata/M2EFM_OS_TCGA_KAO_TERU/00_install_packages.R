# Run once in RStudio with R 4.5.x.

cran_packages <- c(
  "data.table",
  "glmnet",
  "MatrixEQTL",
  "matrixStats",
  "survival"
)

bioc_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "GenomicFeatures",
  "GenomicRanges",
  "IRanges",
  "minfi",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "hugene10sttranscriptcluster.db"
)

missing_cran <- cran_packages[!vapply(
  cran_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

r_minor <- paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
if (identical(r_minor, "4.5")) {
  BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)
}

missing_bioc <- bioc_packages[!vapply(
  bioc_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_bioc) > 0L) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

cat("Package installation/check completed.\n")
cat("R:", R.version.string, "\n")
cat("Bioconductor:", as.character(BiocManager::version()), "\n")
