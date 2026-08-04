#!/usr/bin/env Rscript

# ============================================================
# STEP 10A — COMBAT HARMONIZATION FOR FIGURE 2B
#
# Input:
#   results/09B_external_os_preprocessed.rds
#   results/09B_external_expression_signature_113genes.txt
#
# Output:
#   results/10_external_os_combat_ready.rds
#   results/10_combat_settings.csv
#   results/10_combat_gene_filter_audit.csv
#   results/10_combat_distribution_before_after.csv
#   results/10_combat_batch_mean_distance.csv
#   results/10_combat_pca_coordinates.csv
#   results/10_combat_pca_batch_audit.csv
#   results/10_combat_pca_before_after.png
#   results/10_combat_pca_before_after.pdf
#   results/10_combat_readiness.csv
#   results/10_combat_sessionInfo.txt
#
# Methodological choices:
#   - ComBat is run jointly on TCGA, Terunuma, and Kao.
#   - batch = cohort.
#   - mod = NULL because the paper does not specify a covariate-preserving
#     model matrix for ComBat and cohort is confounded with platform.
#   - par.prior = TRUE, prior.plots = FALSE, mean.only = FALSE.
#   - TDM is not re-applied because the available TCGA matrix is already
#     processed and log-scale; raw RSEM counts are not being reprocessed here.
#
# This script does not train survival models and does not compute rorS.
# rorS will be computed in Step 10B after this audit passes.
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

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_file <- file.path(
  result_dir,
  "09B_external_os_preprocessed.rds"
)

signature_file <- file.path(
  result_dir,
  "09B_external_expression_signature_113genes.txt"
)

required_packages <- c(
  "sva",
  "matrixStats",
  "ggplot2"
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
    ),
    "\nInstall missing packages before running Step 10A."
  )
}

suppressPackageStartupMessages({
  library(sva)
  library(matrixStats)
  library(ggplot2)
})

if (!file.exists(input_file)) {
  stop(
    "Missing input file: ",
    input_file
  )
}

if (!file.exists(signature_file)) {
  stop(
    "Missing signature file: ",
    signature_file
  )
}


# ------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------

clean_text <- function(x) {
  y <- trimws(
    as.character(x)
  )

  y[
    is.na(y) |
      y %in% c(
        "",
        "NA",
        "N/A",
        "NaN",
        "--",
        "."
      )
  ] <- NA_character_

  y
}

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
      " expression lacks row or column names."
    )
  }

  if (
    anyNA(rownames(x)) ||
      anyNA(colnames(x))
  ) {
    stop(
      cohort,
      " expression contains missing row or column names."
    )
  }

  if (anyDuplicated(rownames(x))) {
    stop(
      cohort,
      " expression contains duplicated genes."
    )
  }

  if (anyDuplicated(colnames(x))) {
    stop(
      cohort,
      " expression contains duplicated samples."
    )
  }

  if (
    anyNA(x) ||
      any(!is.finite(x))
  ) {
    stop(
      cohort,
      " expression contains non-finite values."
    )
  }

  x
}

distribution_row <- function(
    x,
    cohort,
    state,
    maximum_values = 1000000L
) {

  total_cells <- as.double(
    length(x)
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

  row_index <- seq.int(
    1L,
    nrow(x),
    by = stride
  )

  column_index <- seq.int(
    1L,
    ncol(x),
    by = stride
  )

  values <- as.numeric(
    x[
      row_index,
      column_index,
      drop = FALSE
    ]
  )

  values <- values[
    is.finite(values)
  ]

  data.frame(
    state = state,
    cohort = cohort,
    genes = nrow(x),
    samples = ncol(x),
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
    sd = stats::sd(values),
    stringsAsFactors = FALSE
  )
}

batch_mean_distance <- function(
    x,
    batch,
    state
) {

  cohort_levels <- levels(batch)

  gene_means <- vapply(
    cohort_levels,
    function(cohort) {
      matrixStats::rowMeans2(
        x[
          ,
          batch == cohort,
          drop = FALSE
        ]
      )
    },
    numeric(
      nrow(x)
    )
  )

  rownames(gene_means) <- rownames(x)

  pairs <- utils::combn(
    cohort_levels,
    2,
    simplify = FALSE
  )

  do.call(
    rbind,
    lapply(
      pairs,
      function(pair) {

        difference <- gene_means[
          ,
          pair[[1]]
        ] -
          gene_means[
            ,
            pair[[2]]
          ]

        data.frame(
          state = state,
          cohort_A = pair[[1]],
          cohort_B = pair[[2]],
          mean_absolute_gene_mean_difference = mean(
            abs(difference)
          ),
          median_absolute_gene_mean_difference = median(
            abs(difference)
          ),
          root_mean_square_gene_mean_difference = sqrt(
            mean(
              difference^2
            )
          ),
          correlation_of_gene_means = suppressWarnings(
            stats::cor(
              gene_means[
                ,
                pair[[1]]
              ],
              gene_means[
                ,
                pair[[2]]
              ],
              method = "pearson"
            )
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}

top_variable_genes <- function(
    x,
    n = 1000L
) {

  variance <- matrixStats::rowVars(x)

  names(variance) <- rownames(x)

  variance <- variance[
    is.finite(variance)
  ]

  ordered <- names(
    sort(
      variance,
      decreasing = TRUE
    )
  )

  ordered[
    seq_len(
      min(
        n,
        length(ordered)
      )
    )
  ]
}

pca_audit <- function(
    x,
    batch,
    state,
    genes,
    components = 10L
) {

  pca_input <- t(
    x[
      genes,
      ,
      drop = FALSE
    ]
  )

  pca <- stats::prcomp(
    pca_input,
    center = TRUE,
    scale. = FALSE,
    rank. = components
  )

  component_count <- min(
    components,
    ncol(pca$x)
  )

  total_variance <- sum(
    apply(
      pca_input,
      2L,
      stats::var
    )
  )

  variance_explained <- (
    pca$sdev[
      seq_len(component_count)
    ]^2
  ) /
    total_variance

  coordinates <- data.frame(
    sample_id = rownames(pca$x),
    cohort = as.character(batch),
    state = state,
    pca$x[
      ,
      seq_len(component_count),
      drop = FALSE
    ],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  audit <- do.call(
    rbind,
    lapply(
      seq_len(component_count),
      function(component) {

        score <- pca$x[
          ,
          component
        ]

        total_ss <- sum(
          (
            score -
              mean(score)
          )^2
        )

        fitted_group_mean <- ave(
          score,
          batch,
          FUN = mean
        )

        between_ss <- sum(
          (
            fitted_group_mean -
              mean(score)
          )^2
        )

        fit <- stats::lm(
          score ~ batch
        )

        anova_table <- stats::anova(fit)

        data.frame(
          state = state,
          component = paste0(
            "PC",
            component
          ),
          variance_explained = variance_explained[
            component
          ],
          batch_R_squared = if (
            total_ss > 0
          ) {
            between_ss / total_ss
          } else {
            NA_real_
          },
          batch_ANOVA_p_value = anova_table[
            1L,
            "Pr(>F)"
          ],
          stringsAsFactors = FALSE
        )
      }
    )
  )

  list(
    coordinates = coordinates,
    audit = audit
  )
}


# ------------------------------------------------------------
# 2. Load Step 09B object
# ------------------------------------------------------------

cat(
  "Loading ",
  input_file,
  "...\n",
  sep = ""
)

input <- readRDS(
  input_file
)

cohorts <- c(
  "TCGA",
  "Terunuma",
  "Kao"
)

required_names <- c(
  cohorts,
  "common_genes"
)

missing_names <- setdiff(
  required_names,
  names(input)
)

if (length(missing_names) > 0L) {
  stop(
    "Input object lacks: ",
    paste(
      missing_names,
      collapse = ", "
    )
  )
}

expression <- lapply(
  cohorts,
  function(cohort) {
    validate_expression(
      input[[cohort]]$expression,
      cohort
    )
  }
)

names(expression) <- cohorts

clinical <- lapply(
  cohorts,
  function(cohort) {
    x <- input[[cohort]]$clinical

    if (is.null(rownames(x))) {
      stop(
        cohort,
        " clinical table lacks row names."
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

common_genes <- unique(
  clean_text(
    input$common_genes
  )
)

common_genes <- common_genes[
  !is.na(common_genes)
]

if (length(common_genes) == 0L) {
  stop(
    "The common gene set is empty."
  )
}

for (cohort in cohorts) {

  absent <- setdiff(
    common_genes,
    rownames(
      expression[[cohort]]
    )
  )

  if (length(absent) > 0L) {
    stop(
      cohort,
      " is missing ",
      length(absent),
      " genes from input$common_genes."
    )
  }

  expression[[cohort]] <- expression[[cohort]][
    common_genes,
    ,
    drop = FALSE
  ]
}

signature_113 <- unique(
  clean_text(
    readLines(
      signature_file,
      warn = FALSE
    )
  )
)

signature_113 <- signature_113[
  !is.na(signature_113)
]

if (length(signature_113) != 113L) {
  stop(
    "Expected 113 external signature genes but found ",
    length(signature_113),
    "."
  )
}

missing_signature <- setdiff(
  signature_113,
  common_genes
)

if (length(missing_signature) > 0L) {
  stop(
    "Signature genes missing from the common gene set: ",
    paste(
      missing_signature,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 3. Combine data and remove zero-variance genes
# ------------------------------------------------------------

sample_counts <- vapply(
  expression,
  ncol,
  integer(1)
)

batch <- factor(
  rep(
    cohorts,
    times = sample_counts
  ),
  levels = cohorts
)

combined_before <- do.call(
  cbind,
  expression
)

all_sample_ids <- colnames(
  combined_before
)

if (anyDuplicated(all_sample_ids)) {
  stop(
    "Sample IDs are duplicated across cohorts."
  )
}

overall_variance <- matrixStats::rowVars(
  combined_before
)

per_cohort_variance <- vapply(
  cohorts,
  function(cohort) {
    matrixStats::rowVars(
      expression[[cohort]]
    )
  },
  numeric(
    nrow(combined_before)
  )
)

rownames(per_cohort_variance) <- rownames(
  combined_before
)

remove_nonfinite <- !is.finite(
  overall_variance
) |
  apply(
    per_cohort_variance,
    1L,
    function(values) {
      any(!is.finite(values))
    }
  )

remove_zero_overall <- overall_variance <= 0

remove_zero_any_cohort <- apply(
  per_cohort_variance,
  1L,
  function(values) {
    any(values <= 0)
  }
)

remove_gene <- (
  remove_nonfinite |
    remove_zero_overall |
    remove_zero_any_cohort
)

removed_genes <- rownames(
  combined_before
)[
  remove_gene
]

removed_signature_genes <- intersect(
  signature_113,
  removed_genes
)

gene_filter_audit <- data.frame(
  metric = c(
    "input_common_genes",
    "nonfinite_variance_genes",
    "zero_variance_overall_genes",
    "zero_variance_in_at_least_one_cohort_genes",
    "total_genes_removed",
    "genes_retained_for_ComBat",
    "signature_genes_before_filter",
    "signature_genes_removed",
    "signature_genes_retained"
  ),
  value = c(
    nrow(combined_before),
    sum(remove_nonfinite),
    sum(remove_zero_overall),
    sum(remove_zero_any_cohort),
    length(removed_genes),
    nrow(combined_before) -
      length(removed_genes),
    length(signature_113),
    length(removed_signature_genes),
    length(
      setdiff(
        signature_113,
        removed_genes
      )
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  gene_filter_audit,
  file.path(
    result_dir,
    "10_combat_gene_filter_audit.csv"
  ),
  row.names = FALSE
)

writeLines(
  removed_genes,
  file.path(
    result_dir,
    "10_combat_removed_genes.txt"
  )
)

writeLines(
  removed_signature_genes,
  file.path(
    result_dir,
    "10_combat_removed_signature_genes.txt"
  )
)

if (length(removed_signature_genes) > 0L) {
  stop(
    "Variance filtering would remove signature genes: ",
    paste(
      removed_signature_genes,
      collapse = ", "
    )
  )
}

if (length(removed_genes) > 0L) {

  retained_genes <- rownames(
    combined_before
  )[
    !remove_gene
  ]

  combined_before <- combined_before[
    retained_genes,
    ,
    drop = FALSE
  ]

  common_genes <- retained_genes

  for (cohort in cohorts) {
    expression[[cohort]] <- expression[[cohort]][
      retained_genes,
      ,
      drop = FALSE
    ]
  }
}


# ------------------------------------------------------------
# 4. Pre-ComBat audits
# ------------------------------------------------------------

distribution_before <- do.call(
  rbind,
  lapply(
    cohorts,
    function(cohort) {
      distribution_row(
        expression[[cohort]],
        cohort,
        "Before_ComBat"
      )
    }
  )
)

distance_before <- batch_mean_distance(
  combined_before,
  batch,
  "Before_ComBat"
)


# ------------------------------------------------------------
# 5. Run ComBat
# ------------------------------------------------------------

combat_settings <- data.frame(
  parameter = c(
    "input_file",
    "genes",
    "samples_total",
    "TCGA_samples",
    "Terunuma_samples",
    "Kao_samples",
    "batch",
    "mod",
    "par.prior",
    "prior.plots",
    "mean.only",
    "TDM_reapplied",
    "note"
  ),
  value = c(
    normalizePath(
      input_file,
      winslash = "/",
      mustWork = TRUE
    ),
    nrow(combined_before),
    ncol(combined_before),
    sample_counts[["TCGA"]],
    sample_counts[["Terunuma"]],
    sample_counts[["Kao"]],
    "cohort",
    "NULL",
    "TRUE",
    "FALSE",
    "FALSE",
    "FALSE",
    paste(
      "ComBat applied jointly to the three available processed",
      "log-scale expression matrices."
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  combat_settings,
  file.path(
    result_dir,
    "10_combat_settings.csv"
  ),
  row.names = FALSE
)

cat(
  "\nRunning ComBat on ",
  nrow(combined_before),
  " genes x ",
  ncol(combined_before),
  " samples...\n",
  sep = ""
)

combat_start <- Sys.time()

combined_after <- sva::ComBat(
  dat = combined_before,
  batch = batch,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE
)

combat_end <- Sys.time()

combined_after <- as.matrix(
  combined_after
)

storage.mode(combined_after) <- "double"

rownames(combined_after) <- rownames(
  combined_before
)

colnames(combined_after) <- colnames(
  combined_before
)

if (
  anyNA(combined_after) ||
    any(!is.finite(combined_after))
) {
  stop(
    "ComBat output contains non-finite values."
  )
}

if (!identical(
  dim(combined_before),
  dim(combined_after)
)) {
  stop(
    "ComBat changed matrix dimensions."
  )
}

combat_minutes <- as.numeric(
  difftime(
    combat_end,
    combat_start,
    units = "mins"
  )
)

cat(
  "ComBat completed in ",
  round(
    combat_minutes,
    3
  ),
  " minutes.\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Split corrected data back into cohorts
# ------------------------------------------------------------

sample_indices <- split(
  seq_len(
    ncol(combined_after)
  ),
  batch
)

expression_after <- lapply(
  cohorts,
  function(cohort) {
    combined_after[
      ,
      sample_indices[[cohort]],
      drop = FALSE
    ]
  }
)

names(expression_after) <- cohorts

for (cohort in cohorts) {

  if (!identical(
    colnames(expression_after[[cohort]]),
    rownames(clinical[[cohort]])
  )) {
    stop(
      "Post-ComBat expression/clinical mismatch for ",
      cohort,
      "."
    )
  }
}


# ------------------------------------------------------------
# 7. Post-ComBat distribution and distance audits
# ------------------------------------------------------------

distribution_after <- do.call(
  rbind,
  lapply(
    cohorts,
    function(cohort) {
      distribution_row(
        expression_after[[cohort]],
        cohort,
        "After_ComBat"
      )
    }
  )
)

distribution_before_after <- rbind(
  distribution_before,
  distribution_after
)

write.csv(
  distribution_before_after,
  file.path(
    result_dir,
    "10_combat_distribution_before_after.csv"
  ),
  row.names = FALSE
)

distance_after <- batch_mean_distance(
  combined_after,
  batch,
  "After_ComBat"
)

distance_table <- rbind(
  distance_before,
  distance_after
)

write.csv(
  distance_table,
  file.path(
    result_dir,
    "10_combat_batch_mean_distance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. PCA audit using the same top-variable genes before/after
# ------------------------------------------------------------

pca_genes <- top_variable_genes(
  combined_before,
  n = 1000L
)

writeLines(
  pca_genes,
  file.path(
    result_dir,
    "10_combat_pca_genes.txt"
  )
)

cat(
  "Running PCA audit on ",
  length(pca_genes),
  " genes...\n",
  sep = ""
)

pca_before <- pca_audit(
  combined_before,
  batch,
  "Before_ComBat",
  pca_genes,
  components = 10L
)

pca_after <- pca_audit(
  combined_after,
  batch,
  "After_ComBat",
  pca_genes,
  components = 10L
)

pca_coordinates <- rbind(
  pca_before$coordinates,
  pca_after$coordinates
)

pca_batch_audit <- rbind(
  pca_before$audit,
  pca_after$audit
)

write.csv(
  pca_coordinates,
  file.path(
    result_dir,
    "10_combat_pca_coordinates.csv"
  ),
  row.names = FALSE
)

write.csv(
  pca_batch_audit,
  file.path(
    result_dir,
    "10_combat_pca_batch_audit.csv"
  ),
  row.names = FALSE
)

plot_data <- pca_coordinates[
  ,
  c(
    "sample_id",
    "cohort",
    "state",
    "PC1",
    "PC2"
  )
]

plot_data$state <- factor(
  plot_data$state,
  levels = c(
    "Before_ComBat",
    "After_ComBat"
  ),
  labels = c(
    "Before ComBat",
    "After ComBat"
  )
)

pca_plot <- ggplot(
  plot_data,
  aes(
    x = PC1,
    y = PC2,
    color = cohort,
    shape = cohort
  )
) +
  geom_point(
    alpha = 0.62,
    size = 1.5
  ) +
  facet_wrap(
    ~ state,
    scales = "free"
  ) +
  labs(
    x = "PC1",
    y = "PC2",
    color = "Cohort",
    shape = "Cohort"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill = "white",
      color = "black"
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    result_dir,
    "10_combat_pca_before_after.png"
  ),
  plot = pca_plot,
  width = 10.5,
  height = 5.2,
  units = "in",
  dpi = 320,
  bg = "white"
)

ggsave(
  file.path(
    result_dir,
    "10_combat_pca_before_after.pdf"
  ),
  plot = pca_plot,
  width = 10.5,
  height = 5.2,
  units = "in",
  device = "pdf",
  bg = "white"
)


# ------------------------------------------------------------
# 9. Save model-ready object for Step 10B and Step 11
# ------------------------------------------------------------

signature_final <- signature_113[
  signature_113 %in%
    rownames(combined_after)
]

if (length(signature_final) != 113L) {
  stop(
    "The final external signature does not contain 113 genes."
  )
}

ready_object <- list(
  settings = list(
    time_column = input$settings$time_column,
    event_column = input$settings$event_column,
    age_column = input$settings$age_column,
    stage_column = input$settings$stage_column,
    ComBat_batch = "cohort",
    ComBat_mod = NULL,
    ComBat_par_prior = TRUE,
    ComBat_prior_plots = FALSE,
    ComBat_mean_only = FALSE,
    TDM_reapplied = FALSE,
    rorS_computed = FALSE,
    NCA_available = FALSE
  ),
  TCGA = list(
    expression = expression_after$TCGA,
    clinical = clinical$TCGA
  ),
  Terunuma = list(
    expression = expression_after$Terunuma,
    clinical = clinical$Terunuma
  ),
  Kao = list(
    expression = expression_after$Kao,
    clinical = clinical$Kao
  ),
  common_genes = rownames(combined_after),
  external_expression_signature = signature_final,
  gene_filter_audit = gene_filter_audit,
  combat_settings = combat_settings,
  distribution_before_after = distribution_before_after,
  batch_mean_distance = distance_table,
  pca_batch_audit = pca_batch_audit,
  combat_elapsed_minutes = combat_minutes,
  created_at = Sys.time()
)

output_rds <- file.path(
  result_dir,
  "10_external_os_combat_ready.rds"
)

saveRDS(
  ready_object,
  output_rds,
  compress = "xz"
)


# ------------------------------------------------------------
# 10. Readiness checks
# ------------------------------------------------------------

mean_distance_before <- mean(
  distance_before$mean_absolute_gene_mean_difference
)

mean_distance_after <- mean(
  distance_after$mean_absolute_gene_mean_difference
)

pc1_before <- pca_batch_audit$batch_R_squared[
  pca_batch_audit$state ==
    "Before_ComBat" &
    pca_batch_audit$component ==
      "PC1"
]

pc1_after <- pca_batch_audit$batch_R_squared[
  pca_batch_audit$state ==
    "After_ComBat" &
    pca_batch_audit$component ==
      "PC1"
]

readiness <- data.frame(
  check = c(
    "Step_09B_object_loaded",
    "identical_gene_order_across_three_cohorts",
    "all_sample_IDs_unique",
    "finite_ComBat_input",
    "ComBat_completed",
    "ComBat_dimensions_preserved",
    "113_gene_signature_preserved",
    "mean_batch_distance_reduced",
    "PC1_batch_R2_reduced",
    "model_ready_object_saved",
    "rorS_computed",
    "NCA_signature_available"
  ),
  status = c(
    TRUE,
    TRUE,
    !anyDuplicated(all_sample_ids),
    all(
      is.finite(combined_before)
    ),
    TRUE,
    identical(
      dim(combined_before),
      dim(combined_after)
    ),
    length(signature_final) == 113L,
    mean_distance_after <
      mean_distance_before,
    length(pc1_before) == 1L &&
      length(pc1_after) == 1L &&
      pc1_after <
        pc1_before,
    file.exists(output_rds),
    FALSE,
    FALSE
  ),
  note = c(
    normalizePath(
      input_file,
      winslash = "/",
      mustWork = TRUE
    ),
    paste0(
      nrow(combined_after),
      " genes."
    ),
    paste0(
      ncol(combined_after),
      " samples."
    ),
    "No missing or infinite values.",
    paste0(
      "Elapsed minutes: ",
      round(
        combat_minutes,
        3
      )
    ),
    paste(
      dim(combined_after),
      collapse = " x "
    ),
    "113 genes retained.",
    paste0(
      signif(
        mean_distance_before,
        5
      ),
      " -> ",
      signif(
        mean_distance_after,
        5
      )
    ),
    paste0(
      signif(
        pc1_before,
        5
      ),
      " -> ",
      signif(
        pc1_after,
        5
      )
    ),
    output_rds,
    "FALSE by design; Step 10B computes rorS.",
    "FALSE; exact NCA signature unavailable."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  readiness,
  file.path(
    result_dir,
    "10_combat_readiness.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Reproducibility report and console output
# ------------------------------------------------------------

end_time <- Sys.time()

total_minutes <- as.numeric(
  difftime(
    end_time,
    start_time,
    units = "mins"
  )
)

capture.output(
  {
    cat("STEP 10A — COMBAT HARMONIZATION\n")
    cat("Created: ", format(end_time), "\n", sep = "")
    cat("Total elapsed minutes: ", total_minutes, "\n\n", sep = "")

    cat("GENE FILTER AUDIT\n")
    print(gene_filter_audit)

    cat("\nCOMBAT SETTINGS\n")
    print(combat_settings)

    cat("\nDISTRIBUTIONS BEFORE/AFTER\n")
    print(distribution_before_after)

    cat("\nBATCH MEAN DISTANCE\n")
    print(distance_table)

    cat("\nPCA BATCH AUDIT\n")
    print(pca_batch_audit)

    cat("\nREADINESS\n")
    print(readiness)

    cat("\nSESSION INFO\n")
    sessionInfo()
  },
  file = file.path(
    result_dir,
    "10_combat_sessionInfo.txt"
  )
)

cat("\n")
cat("============================================================\n")
cat("STEP 10A COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nGENE FILTER AUDIT\n")
print(
  gene_filter_audit,
  row.names = FALSE
)

cat("\nDISTRIBUTIONS BEFORE/AFTER\n")
print(
  distribution_before_after[
    ,
    c(
      "state",
      "cohort",
      "genes",
      "samples",
      "mean",
      "median",
      "sd",
      "q01",
      "q99"
    )
  ],
  row.names = FALSE
)

cat("\nBATCH MEAN DISTANCE\n")
print(
  distance_table[
    ,
    c(
      "state",
      "cohort_A",
      "cohort_B",
      "mean_absolute_gene_mean_difference",
      "correlation_of_gene_means"
    )
  ],
  row.names = FALSE
)

cat("\nPCA BATCH AUDIT — PC1 TO PC3\n")
print(
  pca_batch_audit[
    pca_batch_audit$component %in%
      c(
        "PC1",
        "PC2",
        "PC3"
      ),
    ,
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nREADINESS\n")
print(
  readiness,
  row.names = FALSE
)

cat(
  "\nModel-ready output:\n",
  output_rds,
  "\n",
  sep = ""
)

cat(
  "\nTotal elapsed minutes: ",
  round(
    total_minutes,
    3
  ),
  "\n",
  sep = ""
)

cat("\nFinished successfully.\n")
