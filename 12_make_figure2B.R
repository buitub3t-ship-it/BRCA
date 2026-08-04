#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
result_dir <- file.path(root_dir, "results")
input_file <- file.path(
  result_dir,
  "11_figure2B_split_cindices_100splits.csv"
)

required_packages <- c("data.table", "ggplot2")
ok <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) {
  stop("Missing packages: ", paste(names(ok)[!ok], collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

dataset_order <- c("TCGA", "Terunuma", "Kao")
model_order <- c(
  "Clin",
  "M2EFM Exp+Clin",
  "Cox Exp+Clin",
  "rorS+Clin",
  "M2EFM Exp",
  "Cox Exp",
  "rorS"
)

model_colors <- c(
  "Clin" = "#8C8C8C",
  "M2EFM Exp+Clin" = "#E69F00",
  "Cox Exp+Clin" = "#3BA272",
  "rorS+Clin" = "#E5D431",
  "M2EFM Exp" = "#56B4E9",
  "Cox Exp" = "#49A87A",
  "rorS" = "#C9BE37"
)

x <- fread(input_file, data.table = FALSE, check.names = FALSE)

required_columns <- c("split", "dataset", "model", "c_index")
missing_columns <- setdiff(required_columns, colnames(x))
if (length(missing_columns) > 0L) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

x$split <- as.integer(x$split)
x$dataset <- as.character(x$dataset)
x$model <- as.character(x$model)
x$c_index <- as.numeric(x$c_index)

if (
  anyNA(x) ||
  any(!is.finite(x$c_index)) ||
  any(x$c_index < 0 | x$c_index > 1)
) {
  stop("Invalid or missing input values.")
}

if (!setequal(unique(x$dataset), dataset_order)) {
  stop("Dataset names do not match expected values.")
}

if (!setequal(unique(x$model), model_order)) {
  stop("Model names do not match expected values.")
}

if (anyDuplicated(x[, c("split", "dataset", "model")])) {
  stop("Duplicated split-dataset-model rows.")
}

expected_rows <- 100L * length(dataset_order) * length(model_order)
if (nrow(x) != expected_rows) {
  stop("Expected ", expected_rows, " rows but found ", nrow(x), ".")
}

for (dataset_name in dataset_order) {
  for (model_name in model_order) {
    observed <- sort(
      x$split[
        x$dataset == dataset_name &
          x$model == model_name
      ]
    )
    if (!identical(observed, 1:100)) {
      stop(
        "Expected splits 1:100 for ",
        dataset_name,
        " / ",
        model_name
      )
    }
  }
}

x$dataset <- factor(x$dataset, levels = dataset_order)
x$model <- factor(x$model, levels = model_order)
x <- x[order(x$dataset, x$model, x$split), , drop = FALSE]


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary_rows <- list()
k <- 1L

for (dataset_name in dataset_order) {
  for (model_name in model_order) {
    values <- x$c_index[
      x$dataset == dataset_name &
        x$model == model_name
    ]

    summary_rows[[k]] <- data.frame(
      dataset = dataset_name,
      model = model_name,
      n = length(values),
      mean = mean(values),
      median = median(values),
      sd = sd(values),
      q25 = unname(quantile(values, 0.25)),
      q75 = unname(quantile(values, 0.75)),
      min = min(values),
      max = max(values),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

summary_table <- do.call(rbind, summary_rows)

write.csv(
  summary_table,
  file.path(result_dir, "12_figure2B_summary.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# Paired Wilcoxon tests
# ------------------------------------------------------------

paired_test <- function(dataset_name, model_a, model_b) {
  a <- x[
    x$dataset == dataset_name & x$model == model_a,
    c("split", "c_index"),
    drop = FALSE
  ]
  b <- x[
    x$dataset == dataset_name & x$model == model_b,
    c("split", "c_index"),
    drop = FALSE
  ]

  colnames(a)[2L] <- "A"
  colnames(b)[2L] <- "B"

  paired <- merge(a, b, by = "split", sort = TRUE)

  if (nrow(paired) != 100L || !identical(paired$split, 1:100)) {
    stop("Pairing failed for ", dataset_name, ": ", model_a, " vs ", model_b)
  }

  delta <- paired$A - paired$B

  test <- suppressWarnings(
    wilcox.test(
      paired$A,
      paired$B,
      paired = TRUE,
      alternative = "two.sided",
      exact = FALSE,
      correct = TRUE
    )
  )

  data.frame(
    dataset = dataset_name,
    model_A = model_a,
    model_B = model_b,
    n_pairs = nrow(paired),
    mean_A = mean(paired$A),
    mean_B = mean(paired$B),
    median_A = median(paired$A),
    median_B = median(paired$B),
    mean_delta_A_minus_B = mean(delta),
    median_delta_A_minus_B = median(delta),
    sd_delta = sd(delta),
    proportion_A_greater_B = mean(delta > 0),
    proportion_equal = mean(delta == 0),
    wilcoxon_statistic = unname(test$statistic),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

key_pairs <- list(
  c("M2EFM Exp+Clin", "Clin"),
  c("M2EFM Exp+Clin", "Cox Exp+Clin"),
  c("M2EFM Exp+Clin", "rorS+Clin"),
  c("M2EFM Exp", "Cox Exp"),
  c("M2EFM Exp", "rorS")
)

key_rows <- list()
k <- 1L

for (dataset_name in dataset_order) {
  for (pair in key_pairs) {
    key_rows[[k]] <- paired_test(
      dataset_name,
      pair[[1L]],
      pair[[2L]]
    )
    k <- k + 1L
  }
}

key_tests <- do.call(rbind, key_rows)
key_tests$p_value_BH_global <- p.adjust(key_tests$p_value, method = "BH")
key_tests$p_value_BH_within_dataset <- ave(
  key_tests$p_value,
  key_tests$dataset,
  FUN = function(values) p.adjust(values, method = "BH")
)
key_tests$significant_BH_global_0_05 <-
  key_tests$p_value_BH_global < 0.05

write.csv(
  key_tests,
  file.path(result_dir, "12_figure2B_key_paired_tests.csv"),
  row.names = FALSE
)

all_pairs <- combn(model_order, 2, simplify = FALSE)
all_rows <- list()
k <- 1L

for (dataset_name in dataset_order) {
  for (pair in all_pairs) {
    all_rows[[k]] <- paired_test(
      dataset_name,
      pair[[1L]],
      pair[[2L]]
    )
    k <- k + 1L
  }
}

all_tests <- do.call(rbind, all_rows)
all_tests$p_value_BH_global <- p.adjust(all_tests$p_value, method = "BH")
all_tests$p_value_BH_within_dataset <- ave(
  all_tests$p_value,
  all_tests$dataset,
  FUN = function(values) p.adjust(values, method = "BH")
)
all_tests$significant_BH_global_0_05 <-
  all_tests$p_value_BH_global < 0.05

write.csv(
  all_tests,
  file.path(result_dir, "12_figure2B_all_pairwise_tests.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------

plot_data <- x
plot_data$dataset <- factor(plot_data$dataset, levels = dataset_order)
plot_data$model <- factor(plot_data$model, levels = model_order)

figure <- ggplot(
  plot_data,
  aes(x = model, y = c_index, fill = model)
) +
  geom_hline(
    yintercept = 0.50,
    linewidth = 0.5,
    linetype = "dashed",
    color = "#555555"
  ) +
  geom_hline(
    yintercept = 0.70,
    linewidth = 0.55,
    color = "#D62728"
  ) +
  geom_boxplot(
    width = 0.72,
    linewidth = 0.55,
    outlier.shape = 16,
    outlier.size = 1.2,
    outlier.alpha = 0.70
  ) +
  facet_grid(. ~ dataset) +
  scale_fill_manual(
    values = model_colors,
    breaks = model_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0.45, 0.90),
    breaks = seq(0.5, 0.9, by = 0.1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = NULL,
    y = "C-index",
    fill = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.65),
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.65
    ),
    strip.text = element_text(face = "bold", size = 12),
    axis.title.y = element_text(
      face = "bold",
      margin = margin(r = 9)
    ),
    axis.text.x = element_text(
      angle = 48,
      hjust = 1,
      vjust = 1,
      size = 8.5
    ),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(t = 12, r = 16, b = 10, l = 12)
  )

png_file <- file.path(result_dir, "12_figure2B_7models.png")
pdf_file <- file.path(result_dir, "12_figure2B_7models.pdf")

ggsave(
  png_file,
  plot = figure,
  width = 15.5,
  height = 6.8,
  units = "in",
  dpi = 320,
  bg = "white"
)

ggsave(
  pdf_file,
  plot = figure,
  width = 15.5,
  height = 6.8,
  units = "in",
  device = "pdf",
  bg = "white"
)


# ------------------------------------------------------------
# Audit and report
# ------------------------------------------------------------

plot_audit <- data.frame(
  metric = c(
    "input_file",
    "input_rows",
    "datasets",
    "models_available",
    "splits_per_dataset_model",
    "expected_C_index_values",
    "observed_C_index_values",
    "NCA_models_included",
    "key_paired_tests",
    "all_pairwise_tests",
    "statistical_test",
    "p_value_adjustment"
  ),
  value = c(
    normalizePath(input_file, winslash = "/", mustWork = TRUE),
    nrow(x),
    length(dataset_order),
    length(model_order),
    100,
    expected_rows,
    nrow(x),
    FALSE,
    nrow(key_tests),
    nrow(all_tests),
    "two-sided paired Wilcoxon signed-rank",
    "Benjamini-Hochberg"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  plot_audit,
  file.path(result_dir, "12_figure2B_plot_audit.csv"),
  row.names = FALSE
)

capture.output(
  {
    cat("STEP 12 — FIGURE 2B\n\n")
    cat("SUMMARY\n")
    print(summary_table)
    cat("\nKEY PAIRED TESTS\n")
    print(key_tests)
    cat("\nPLOT AUDIT\n")
    print(plot_audit)
    cat(
      "\nNOTE\nSeven models are plotted. NCA Exp and NCA Exp+Clin ",
      "were not reconstructed because the exact NCA signature ",
      "was unavailable.\n",
      sep = ""
    )
    cat("\nSESSION INFO\n")
    sessionInfo()
  },
  file = file.path(result_dir, "12_figure2B_sessionInfo.txt")
)

cat("\n============================================================\n")
cat("STEP 12 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nFIGURE 2B SUMMARY\n")
print(
  summary_table[, c("dataset", "model", "n", "mean", "median", "sd")],
  row.names = FALSE
)

cat("\nKEY PAIRED TESTS\n")
print(
  key_tests[
    ,
    c(
      "dataset",
      "model_A",
      "model_B",
      "median_delta_A_minus_B",
      "proportion_A_greater_B",
      "p_value",
      "p_value_BH_global"
    )
  ],
  row.names = FALSE
)

cat("\nFigure saved:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat("\nFinished successfully.\n")
