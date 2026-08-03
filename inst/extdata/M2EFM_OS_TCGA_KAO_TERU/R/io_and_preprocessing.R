assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Run 00_install_packages.R first."
    )
  }
}

read_feature_matrix <- function(path, id_col) {
  dt <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  if (!id_col %in% names(dt)) {
    stop("Column '", id_col, "' not found in ", basename(path))
  }

  ids <- as.character(dt[[id_col]])
  if (anyNA(ids) || any(!nzchar(ids))) {
    stop("Missing/blank feature IDs in ", basename(path))
  }
  if (anyDuplicated(ids)) {
    stop("Duplicated feature IDs in ", basename(path))
  }

  dt[[id_col]] <- NULL
  mat <- data.matrix(dt)
  rownames(mat) <- ids
  storage.mode(mat) <- "double"

  if (any(!is.finite(mat))) {
    stop("Non-finite values found in ", basename(path))
  }
  if (anyDuplicated(colnames(mat))) {
    stop("Duplicated sample IDs in ", basename(path))
  }
  mat
}

read_clinical <- function(path) {
  clin <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c("NA", "N/A", "--", "")
  )
  if (!"sample_id" %in% names(clin)) {
    stop("sample_id not found in ", basename(path))
  }
  clin$sample_id <- as.character(clin$sample_id)
  if (anyDuplicated(clin$sample_id)) {
    stop("Duplicated sample_id in ", basename(path))
  }
  rownames(clin) <- clin$sample_id
  clin
}

reduce_stage <- function(x) {
  s <- trimws(as.character(x))
  s[s %in% c("", "NA", "N/A", "Stage X", "X")] <- NA_character_
  s <- sub("^stage\\s*", "", s, ignore.case = TRUE)
  s <- toupper(s)

  out <- rep(NA_character_, length(s))
  out[grepl("^I($|A$|B$)", s)] <- "Stage I"
  out[grepl("^II($|A$|B$)", s)] <- "Stage II"
  out[grepl("^III($|A$|B$|C$)", s)] <- "Stage III"
  out[grepl("^IV($|A$|B$)", s)] <- "Stage IV"
  factor(out, levels = STAGE_LEVELS)
}

clean_clinical <- function(clin, cohort) {
  required <- c(
    "sample_id", "age.Dx", "pathologic_stage",
    "OVERALL.SURVIVAL", "overall.survival.indicator"
  )
  absent <- setdiff(required, names(clin))
  if (length(absent) > 0L) {
    stop(cohort, " clinical file lacks: ", paste(absent, collapse = ", "))
  }

  clin$age.Dx <- suppressWarnings(as.numeric(clin$age.Dx))
  clin$pathologic_stage <- reduce_stage(clin$pathologic_stage)
  clin$overall.survival.indicator <- suppressWarnings(
    as.integer(clin$overall.survival.indicator)
  )
  raw_time <- suppressWarnings(as.numeric(clin$OVERALL.SURVIVAL))

  # Standardize only for reporting/KM/calibration. C-index is invariant to
  # positive unit conversions.
  cohort_upper <- toupper(cohort)
  if (cohort_upper == "TCGA") {
    clin$OS_days <- raw_time
  } else if (cohort_upper == "KAO") {
    clin$OS_days <- raw_time * 365.25
  } else if (cohort_upper %in% c("TERU", "TERUNUMA")) {
    clin$OS_days <- raw_time * 365.25 / 12
  } else {
    stop("Unknown cohort: ", cohort)
  }

  invalid_event <- !is.na(clin$overall.survival.indicator) &
    !clin$overall.survival.indicator %in% c(0L, 1L)
  if (any(invalid_event)) {
    stop(cohort, " has survival indicators outside {0,1}.")
  }

  clin
}

map_terunuma_to_gene <- function(probe_matrix) {
  assert_packages(c("AnnotationDbi", "hugene10sttranscriptcluster.db"))
  chip_db <- get(
    "hugene10sttranscriptcluster.db",
    envir = asNamespace("hugene10sttranscriptcluster.db")
  )

  probe_ids <- rownames(probe_matrix)
  mapping <- suppressMessages(AnnotationDbi::select(
    chip_db,
    keys = probe_ids,
    keytype = "PROBEID",
    columns = c("SYMBOL")
  ))

  mapping$PROBEID <- as.character(mapping$PROBEID)
  mapping$SYMBOL <- trimws(as.character(mapping$SYMBOL))
  mapping <- unique(mapping[
    !is.na(mapping$SYMBOL) & nzchar(mapping$SYMBOL),
    c("PROBEID", "SYMBOL"),
    drop = FALSE
  ])

  # Drop transcript clusters mapping to more than one symbol.
  n_symbol_per_probe <- table(mapping$PROBEID)
  mapping <- mapping[
    n_symbol_per_probe[mapping$PROBEID] == 1L,
    ,
    drop = FALSE
  ]
  symbol_by_probe <- setNames(mapping$SYMBOL, mapping$PROBEID)

  keep_probes <- intersect(probe_ids, names(symbol_by_probe))
  if (length(keep_probes) == 0L) {
    stop("No Terunuma transcript-cluster IDs mapped to gene symbols.")
  }

  row_indices <- match(keep_probes, probe_ids)
  groups <- split(row_indices, symbol_by_probe[keep_probes])

  # Author's preprocessing aggregated probe sets mapping to the same gene
  # using the median expression value.
  gene_matrix <- t(vapply(
    groups,
    function(idx) {
      if (length(idx) == 1L) {
        as.numeric(probe_matrix[idx, ])
      } else {
        apply(probe_matrix[idx, , drop = FALSE], 2L, stats::median, na.rm = TRUE)
      }
    },
    FUN.VALUE = numeric(ncol(probe_matrix))
  ))

  rownames(gene_matrix) <- names(groups)
  colnames(gene_matrix) <- colnames(probe_matrix)
  storage.mode(gene_matrix) <- "double"

  list(
    expression = gene_matrix,
    mapping = mapping,
    n_input_probes = nrow(probe_matrix),
    n_mapped_unambiguous_probes = length(keep_probes),
    n_gene_symbols = nrow(gene_matrix)
  )
}

align_clinical_to_expression <- function(clin, expression) {
  ids <- intersect(colnames(expression), rownames(clin))
  expression <- expression[, ids, drop = FALSE]
  clin <- clin[ids, , drop = FALSE]
  list(expression = expression, clinical = clin)
}

complete_os_rows <- function(clin) {
  complete.cases(
    clin[, c(
      "OS_days", "overall.survival.indicator",
      "age.Dx", "pathologic_stage"
    ), drop = FALSE]
  ) & clin$OS_days >= 0
}

make_audit_row <- function(dataset, matrix, clinical = NULL) {
  data.frame(
    dataset = dataset,
    n_features = nrow(matrix),
    n_molecular_samples = ncol(matrix),
    n_clinical_samples = if (is.null(clinical)) NA_integer_ else nrow(clinical),
    n_id_overlap = if (is.null(clinical)) NA_integer_ else
      length(intersect(colnames(matrix), rownames(clinical))),
    n_complete_os = if (is.null(clinical)) NA_integer_ else
      sum(complete_os_rows(clinical)),
    molecular_missing_cells = sum(is.na(matrix)),
    stringsAsFactors = FALSE
  )
}
