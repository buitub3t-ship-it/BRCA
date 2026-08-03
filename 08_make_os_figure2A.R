#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
OUT <- file.path(ROOT, "results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
INPUT <- file.path(OUT, "07_nine_model_split_cindices_100splits.csv")

pkgs <- c("data.table", "ggplot2")
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) {
  stop(
    "Missing packages: ", paste(names(ok)[!ok], collapse = ", "),
    "\nInstall with: install.packages(c('data.table','ggplot2'), repos='https://cloud.r-project.org')"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!file.exists(INPUT)) stop("Missing input: ", INPUT)

MODEL_ORDER <- c(
  "Clin",
  "M2EFM Meth+Exp+Clin",
  "M2EFM Exp+Clin",
  "Cox Meth+Exp+Clin",
  "rorS+Clin",
  "M2EFM Meth+Exp",
  "M2EFM Exp",
  "Cox Meth+Exp",
  "rorS"
)

MODEL_COLORS <- c(
  "Clin" = "#8C8C8C",
  "M2EFM Meth+Exp+Clin" = "#E69F00",
  "M2EFM Exp+Clin" = "#56B4E9",
  "Cox Meth+Exp+Clin" = "#3BA272",
  "rorS+Clin" = "#E5D431",
  "M2EFM Meth+Exp" = "#D89000",
  "M2EFM Exp" = "#72B7E2",
  "Cox Meth+Exp" = "#49A87A",
  "rorS" = "#D8CD3F"
)

PUBLISHED_MEDIANS <- c(
  "Clin" = 0.753,
  "M2EFM Meth+Exp+Clin" = 0.790,
  "M2EFM Exp+Clin" = 0.787,
  "Cox Meth+Exp+Clin" = 0.728,
  "rorS+Clin" = 0.775,
  "M2EFM Meth+Exp" = 0.688,
  "M2EFM Exp" = 0.682,
  "Cox Meth+Exp" = 0.662,
  "rorS" = 0.636
)

x <- fread(INPUT, data.table = FALSE, check.names = FALSE)
need <- c("split", "model", "c_index")
miss <- setdiff(need, names(x))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

x$split <- suppressWarnings(as.integer(x$split))
x$model <- as.character(x$model)
x$c_index <- suppressWarnings(as.numeric(x$c_index))

if (anyNA(x[, need])) stop("NA/non-convertible values detected.")
if (any(!is.finite(x$c_index)) || any(x$c_index < 0 | x$c_index > 1)) {
  stop("C-index must be finite and in [0,1].")
}
if (length(setdiff(unique(x$model), MODEL_ORDER))) {
  stop("Unexpected model names: ", paste(setdiff(unique(x$model), MODEL_ORDER), collapse = ", "))
}
if (length(setdiff(MODEL_ORDER, unique(x$model)))) {
  stop("Missing models: ", paste(setdiff(MODEL_ORDER, unique(x$model)), collapse = ", "))
}
if (anyDuplicated(x[, c("split", "model")])) stop("Duplicate split-model rows detected.")
if (nrow(x) != 900L) stop("Expected 900 rows; found ", nrow(x))

counts <- table(factor(x$model, levels = MODEL_ORDER))
if (any(counts != 100L)) {
  stop("Each model must have 100 values: ", paste(names(counts), counts, sep = "=", collapse = ", "))
}
for (m in MODEL_ORDER) {
  if (!identical(sort(x$split[x$model == m]), 1:100)) {
    stop("Model does not contain exactly splits 1:100: ", m)
  }
}

x$model <- factor(x$model, levels = MODEL_ORDER)
x <- x[order(x$model, x$split), , drop = FALSE]

summarize_model <- function(v, m) {
  data.frame(
    model = m,
    n = length(v),
    mean = mean(v),
    median = median(v),
    sd = sd(v),
    q25 = unname(quantile(v, 0.25)),
    q75 = unname(quantile(v, 0.75)),
    min = min(v),
    max = max(v),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, lapply(MODEL_ORDER, function(m) {
  summarize_model(x$c_index[x$model == m], m)
}))
write.csv(summary_table, file.path(OUT, "08_os_nine_model_summary.csv"), row.names = FALSE)

summary_vs_published <- summary_table
summary_vs_published$published_median <- as.numeric(PUBLISHED_MEDIANS[summary_vs_published$model])
summary_vs_published$replication_minus_published_median <-
  summary_vs_published$median - summary_vs_published$published_median
write.csv(
  summary_vs_published,
  file.path(OUT, "08_os_summary_vs_published.csv"),
  row.names = FALSE
)

wide <- reshape(
  x[, c("split", "model", "c_index")],
  idvar = "split",
  timevar = "model",
  direction = "wide"
)
names(wide) <- sub("^c_index\\.", "", names(wide))
wide <- wide[order(wide$split), , drop = FALSE]
if (!identical(wide$split, 1:100)) stop("Paired table lost split order.")

paired_test <- function(a_name, b_name, alternative = "two.sided") {
  a <- as.numeric(wide[[a_name]])
  b <- as.numeric(wide[[b_name]])
  d <- a - b
  wt <- suppressWarnings(wilcox.test(
    a, b, paired = TRUE, alternative = alternative,
    exact = FALSE, correct = TRUE
  ))
  data.frame(
    model_A = a_name,
    model_B = b_name,
    alternative = alternative,
    n_pairs = length(d),
    mean_A = mean(a),
    mean_B = mean(b),
    median_A = median(a),
    median_B = median(b),
    mean_delta_A_minus_B = mean(d),
    median_delta_A_minus_B = median(d),
    sd_delta = sd(d),
    q025_delta = unname(quantile(d, 0.025)),
    q975_delta = unname(quantile(d, 0.975)),
    proportion_A_greater_B = mean(d > 0),
    proportion_equal = mean(d == 0),
    wilcoxon_statistic = unname(wt$statistic),
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
}

key_pairs <- list(
  c("M2EFM Meth+Exp+Clin", "rorS+Clin"),
  c("M2EFM Meth+Exp+Clin", "Clin"),
  c("M2EFM Meth+Exp+Clin", "Cox Meth+Exp+Clin"),
  c("M2EFM Meth+Exp+Clin", "M2EFM Exp+Clin"),
  c("M2EFM Meth+Exp", "rorS"),
  c("M2EFM Meth+Exp", "Cox Meth+Exp"),
  c("M2EFM Meth+Exp", "M2EFM Exp")
)
key_tests <- do.call(rbind, lapply(key_pairs, function(z) {
  paired_test(z[1], z[2], alternative = "two.sided")
}))
key_tests$p_value_BH <- p.adjust(key_tests$p_value, method = "BH")
key_tests$significant_BH_0_05 <- key_tests$p_value_BH < 0.05
write.csv(key_tests, file.path(OUT, "08_os_key_paired_contrasts.csv"), row.names = FALSE)

all_pairs <- combn(MODEL_ORDER, 2, simplify = FALSE)
all_tests <- do.call(rbind, lapply(all_pairs, function(z) {
  paired_test(z[1], z[2], alternative = "two.sided")
}))
all_tests$p_value_BH <- p.adjust(all_tests$p_value, method = "BH")
all_tests$significant_BH_0_05 <- all_tests$p_value_BH < 0.05
write.csv(all_tests, file.path(OUT, "08_os_all_pairwise_paired_tests.csv"), row.names = FALSE)

p <- ggplot(x, aes(x = model, y = c_index, fill = model)) +
  geom_hline(yintercept = 0.70, linewidth = 0.65, color = "#D62728") +
  geom_boxplot(
    width = 0.72, linewidth = 0.55,
    outlier.shape = 16, outlier.size = 1.45, outlier.alpha = 0.75
  ) +
  scale_fill_manual(values = MODEL_COLORS, breaks = MODEL_ORDER, drop = FALSE) +
  scale_y_continuous(
    limits = c(0.45, 0.93),
    breaks = seq(0.5, 0.9, by = 0.1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(x = NULL, y = "C-index", fill = "Model") +
  annotate("text", x = 0.55, y = 0.925, label = "A", fontface = "bold", size = 5, hjust = 0) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linewidth = 0.65, color = "black"),
    axis.title.y = element_text(face = "bold", margin = margin(r = 9)),
    axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8.5),
    legend.key.height = grid::unit(0.42, "cm"),
    plot.margin = margin(t = 12, r = 16, b = 10, l = 12)
  )

png_file <- file.path(OUT, "08_os_figure2A_replication.png")
pdf_file <- file.path(OUT, "08_os_figure2A_replication.pdf")

ggsave(png_file, p, width = 13.5, height = 6.8, units = "in", dpi = 320, bg = "white")
ggsave(pdf_file, p, width = 13.5, height = 6.8, units = "in", device = "pdf", bg = "white")

audit <- data.frame(
  metric = c(
    "input_rows", "models", "splits_per_model", "total_C_indices",
    "minimum_C_index", "maximum_C_index", "reference_line_C_index",
    "key_two_sided_tests", "all_two_sided_pairwise_tests"
  ),
  value = c(
    nrow(x), length(MODEL_ORDER), 100, nrow(x), min(x$c_index), max(x$c_index),
    0.70, nrow(key_tests), nrow(all_tests)
  ),
  stringsAsFactors = FALSE
)
write.csv(audit, file.path(OUT, "08_os_plot_audit.csv"), row.names = FALSE)

capture.output({
  cat("TCGA internal OS Figure 2A replication\n\n")
  cat("Input: ", INPUT, "\n", sep = "")
  cat("Plot uses this replication's 900 C-index values.\n")
  cat("Tests are paired Wilcoxon tests across the same 100 splits.\n")
  cat("They are not the paper's random-gene-set bootstrap p-values.\n\n")
  print(summary_table)
  cat("\nKey paired contrasts:\n")
  print(key_tests)
  cat("\nAudit:\n")
  print(audit)
  cat("\n")
  sessionInfo()
}, file = file.path(OUT, "08_os_figure2A_sessionInfo.txt"))

cat("\nNINE-MODEL OS SUMMARY\n")
print(summary_table[, c("model", "n", "mean", "median", "sd")], row.names = FALSE)

cat("\nKEY TWO-SIDED PAIRED CONTRASTS\n")
print(
  key_tests[, c(
    "model_A", "model_B", "median_delta_A_minus_B",
    "proportion_A_greater_B", "p_value", "p_value_BH"
  )],
  row.names = FALSE
)

cat("\nFigure saved:\n", png_file, "\n", pdf_file, "\n", sep = "")
cat("\nFinished successfully.\n")
