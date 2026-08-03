beta_to_m <- function(beta, eps = 1e-6) {
  beta <- pmin(pmax(beta, eps), 1 - eps)
  log2(beta / (1 - beta))
}

get_probe_positions_hg19 <- function(probes) {
  assert_packages(c(
    "minfi",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  ))
  
  pkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  
  # Annotation object hiện cần package được attach vào search list.
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
  
  annotation_object <- get(
    pkg,
    envir = as.environment(paste0("package:", pkg))
  )
  
  annotation <- as.data.frame(
    minfi::getAnnotation(annotation_object)
  )
  
  idx <- match(probes, rownames(annotation))
  
  out <- data.frame(
    snpid = probes,
    chr = as.character(annotation$chr[idx]),
    pos = suppressWarnings(as.integer(annotation$pos[idx])),
    stringsAsFactors = FALSE
  )
  
  out <- out[complete.cases(out), , drop = FALSE]
  out <- out[
    grepl("^chr([0-9]+|X|Y)$", out$chr),
    ,
    drop = FALSE
  ]
  
  out
}

map_symbol_to_entrez_with_alias <- function(symbols) {
  assert_packages(c("AnnotationDbi", "org.Hs.eg.db"))
  org_db <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))

  direct <- suppressMessages(AnnotationDbi::mapIds(
    org_db,
    keys = symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  ))

  missing <- names(direct)[is.na(direct)]
  if (length(missing) > 0L) {
    alias_map <- suppressMessages(AnnotationDbi::select(
      org_db,
      keys = missing,
      columns = "ENTREZID",
      keytype = "ALIAS"
    ))
    alias_map <- alias_map[!is.na(alias_map$ENTREZID), , drop = FALSE]
    alias_map <- alias_map[!duplicated(alias_map$ALIAS), , drop = FALSE]
    direct[alias_map$ALIAS] <- alias_map$ENTREZID
  }
  direct
}

get_gene_positions_hg19 <- function(symbols) {
  assert_packages(c(
    "GenomicFeatures", "GenomicRanges", "IRanges",
    "TxDb.Hsapiens.UCSC.hg19.knownGene"
  ))

  entrez <- map_symbol_to_entrez_with_alias(symbols)
  txdb <- get(
    "TxDb.Hsapiens.UCSC.hg19.knownGene",
    envir = asNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")
  )
  gr <- GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)
  gr_df <- data.frame(
    ENTREZID = names(gr),
    chr = as.character(GenomicRanges::seqnames(gr)),
    left = as.integer(IRanges::start(gr)),
    right = as.integer(IRanges::end(gr)),
    stringsAsFactors = FALSE
  )
  rownames(gr_df) <- gr_df$ENTREZID

  idx <- match(as.character(entrez), gr_df$ENTREZID)
  out <- data.frame(
    geneid = symbols,
    chr = gr_df$chr[idx],
    left = gr_df$left[idx],
    right = gr_df$right[idx],
    stringsAsFactors = FALSE
  )
  out <- out[complete.cases(out), , drop = FALSE]
  out <- out[grepl("^chr([0-9]+|X|Y)$", out$chr), , drop = FALSE]
  out <- out[!duplicated(out$geneid), , drop = FALSE]
  out
}

new_sliced_data <- function(mat) {
  sliced_generator <- get("SlicedData", envir = asNamespace("MatrixEQTL"))
  out <- sliced_generator$new()
  out$CreateFromMatrix(as.matrix(mat))
  out
}

run_m2eqtl_discovery <- function(
    tcga_expression,
    tcga_methylation_beta,
    probelist,
    num_probes = NUM_DISCOVERY_PROBES,
    num_trans = NUM_TRANS_GENES,
    cis_p = CIS_P_THRESHOLD,
    trans_p = TRANS_P_THRESHOLD,
    cis_distance = CIS_DISTANCE_BP,
    mad_threshold = MAD_THRESHOLD) {

  assert_packages(c("MatrixEQTL", "matrixStats"))

  tumor_exp_ids <- grep("-01$", colnames(tcga_expression), value = TRUE)
  tumor_meth_ids <- grep("-01$", colnames(tcga_methylation_beta), value = TRUE)
  common_samples <- intersect(tumor_exp_ids, tumor_meth_ids)
  if (length(common_samples) < 100L) {
    stop("Too few matched TCGA tumor expression/methylation samples.")
  }

  selected_probes <- probelist[probelist %in% rownames(tcga_methylation_beta)]
  selected_probes <- head(selected_probes, num_probes)
  if (length(selected_probes) < num_probes) {
    stop("Only ", length(selected_probes), " of the requested ", num_probes,
         " probes are available.")
  }

  exp_mat <- tcga_expression[, common_samples, drop = FALSE]
  meth_beta <- tcga_methylation_beta[selected_probes, common_samples, drop = FALSE]

  gene_mad <- matrixStats::rowMads(exp_mat)
  exp_mat <- exp_mat[is.finite(gene_mad) & gene_mad > mad_threshold, , drop = FALSE]

  probe_pos <- get_probe_positions_hg19(rownames(meth_beta))
  gene_pos <- get_gene_positions_hg19(rownames(exp_mat))

  meth_beta <- meth_beta[probe_pos$snpid, , drop = FALSE]
  exp_mat <- exp_mat[gene_pos$geneid, , drop = FALSE]

  if (!identical(colnames(meth_beta), colnames(exp_mat))) {
    stop("Expression and methylation sample ordering failed.")
  }

  snps <- new_sliced_data(beta_to_m(meth_beta))
  genes <- new_sliced_data(exp_mat)
  sliced_generator <- get("SlicedData", envir = asNamespace("MatrixEQTL"))
  covariates <- sliced_generator$new()

  output_trans <- tempfile(fileext = ".txt")
  output_cis <- tempfile(fileext = ".txt")
  on.exit(unlink(c(output_trans, output_cis)), add = TRUE)

  me <- MatrixEQTL::Matrix_eQTL_main(
    snps = snps,
    gene = genes,
    cvrt = covariates,
    output_file_name = output_trans,
    pvOutputThreshold = trans_p,
    useModel = get("modelLINEAR", envir = asNamespace("MatrixEQTL")),
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = output_cis,
    pvOutputThreshold.cis = cis_p,
    snpspos = probe_pos,
    genepos = gene_pos,
    cisDist = cis_distance,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )

  cis <- as.data.frame(me$cis$eqtls)
  trans <- as.data.frame(me$trans$eqtls)
  if (nrow(trans) == 0L) {
    stop("No trans m2eQTL passed p < ", trans_p, ".")
  }

  cis <- cis[order(-abs(cis$beta)), , drop = FALSE]
  trans <- trans[order(-abs(trans$beta)), , drop = FALSE]

  # Author code: one strongest cis association per probe; one strongest trans
  # association per gene, then keep the top num_trans genes by |beta|.
  cis_unique_probe <- cis[!duplicated(cis$snps), , drop = FALSE]
  trans_unique_gene <- trans[!duplicated(trans$gene), , drop = FALSE]
  trans_selected <- head(trans_unique_gene, num_trans)

  selected_m2e_probes <- unique(as.character(trans_selected$snps))
  cis_replacements <- cis_unique_probe[
    cis_unique_probe$snps %in% selected_m2e_probes,
    ,
    drop = FALSE
  ]

  trans_genes <- unique(as.character(trans_selected$gene))
  cis_genes <- unique(as.character(cis_replacements$gene))
  expression_signature <- unique(c(trans_genes, cis_genes))

  list(
    matrixeqtl = me,
    cis_all_significant = cis,
    trans_all_significant = trans,
    cis_unique_probe = cis_unique_probe,
    trans_selected = trans_selected,
    cis_replacements = cis_replacements,
    trans_genes = trans_genes,
    cis_genes = cis_genes,
    expression_signature = expression_signature,
    methylation_probes = selected_m2e_probes,
    matched_samples = common_samples,
    discovery_probes = selected_probes,
    probe_positions = probe_pos,
    gene_positions = gene_pos
  )
}
