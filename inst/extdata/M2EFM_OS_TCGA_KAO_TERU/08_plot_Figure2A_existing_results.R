# =============================================================================
# 08_plot_Figure2A_existing_results.R
#
# Doc ket qua THUC TE da duoc tao boi 03_train_internal_full.R.
# Khong su dung file khong ton tai: results/03_internal_Cindex_9models.csv
#
# Project:
# D:/M2EFM_BRCA_Replication/M2EFM-master/inst/extdata/M2EFM_OS_TCGA_KAO_TERU
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

required <- c("data.table", "ggplot2")
missing_pkg <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkg) > 0L) {
  stop(
    "Thieu package: ", paste(missing_pkg, collapse = ", "),
    "\nChay: install.packages(c(\"data.table\", \"ggplot2\"))"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. PATH DUNG THEO THU MUC TRONG ANH
# -----------------------------------------------------------------------------

PROJECT_DIR <- "D:/M2EFM_BRCA_Replication/M2EFM-master/inst/extdata/M2EFM_OS_TCGA_KAO_TERU"
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)
RESULT_DIR  <- file.path(PROJECT_DIR, "results")

if (!dir.exists(RESULT_DIR)) {
  stop("Khong tim thay thu muc results: ", RESULT_DIR)
}

setwd(PROJECT_DIR)

# -----------------------------------------------------------------------------
# 2. TIM FILE C-INDEX THUC TE
# Uu tien file 100 random splits dang hien trong anh.
# -----------------------------------------------------------------------------

preferred_names <- c(
  "03_internal_full_monte_carlo_metrics.csv",
  "03_internal_full_metrics.csv",
  "03_final_internal_full_metrics.csv"
)

all_csv <- list.files(
  RESULT_DIR,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

# Khong doc output cua script nay neu chay lai.
all_csv <- all_csv[!grepl("^08_", basename(all_csv))]

preferred_paths <- file.path(RESULT_DIR, preferred_names)
preferred_paths <- preferred_paths[file.exists(preferred_paths)]

# Neu ten file tren Windows bi an phan .csv hoac co khac nhe, tim theo pattern.
pattern_hits <- all_csv[
  grepl("03.*internal.*(monte|metric)", basename(all_csv), ignore.case = TRUE)
]

candidate_files <- unique(c(preferred_paths, pattern_hits))

if (length(candidate_files) == 0L) {
  stop(
    "Khong tim thay file metrics cua buoc 03 trong:\n", RESULT_DIR,
    "\n\nCac CSV dang co:\n- ",
    paste(basename(all_csv), collapse = "\n- ")
  )
}

# Chon file co nhieu dong nhat, vi file Monte Carlo thuong co 100 splits x models.
row_counts <- vapply(candidate_files, function(f) {
  x <- tryCatch(fread(f, nrows = Inf, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x)) return(-1L)
  nrow(x)
}, integer(1))

INPUT_FILE <- candidate_files[which.max(row_counts)]
cat("Input file :", INPUT_FILE, "\n")
cat("Output dir:", RESULT_DIR, "\n\n")

raw <- fread(
  INPUT_FILE,
  data.table = FALSE,
  na.strings = c("NA", "N/A", "--", "")
)

cat("Columns found:\n")
print(names(raw))
cat("\nFirst rows:\n")
print(utils::head(raw))

# -----------------------------------------------------------------------------
# 3. HAM NHAN DIEN COT VA TEN MODEL
# -----------------------------------------------------------------------------

key <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))

find_col <- function(nms, candidates) {
  idx <- match(key(candidates), key(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) return(NA_character_)
  nms[idx[1L]]
}

canonical_model <- function(x) {
  k <- key(x)
  out <- as.character(x)

  out[k %in% c("clin", "clinical", "clinicalonly")] <- "Clin"

  # Trong file internal_full, M2EFM_ExpClin la model tich hop Meth+Exp+Clin.
  out[k %in% c(
    "m2efmexpclin", "m2efmmethexpclin", "m2efmintegratedclin",
    "m2efmmethplusexpplusclin"
  )] <- "M2EFM Meth+Exp+Clin"

  out[k %in% c(
    "m2efmmolecular", "m2efmmethexp", "m2efmintegrated",
    "m2efmmethplusexp"
  )] <- "M2EFM Meth+Exp"

  out[k %in% c("m2efmexpressionclin", "m2efmexpressionplusclin")] <- "M2EFM Exp+Clin"
  out[k %in% c("m2efmexpression", "m2efmexp") ] <- "M2EFM Exp"

  out[k %in% c("coxmethexpclin", "coxridgemethexpclin")] <- "Cox Meth+Exp+Clin"
  out[k %in% c("coxmethexp", "coxridgemethexp")] <- "Cox Meth+Exp"

  out[k %in% c("rorsclin", "pam50rorsclin")] <- "rorS+Clin"
  out[k %in% c("rors", "pam50rors")] <- "rorS"

  out
}

MODEL_LEVELS <- c(
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

model_col <- find_col(names(raw), c("model", "method", "model_name", "approach"))
cindex_col <- find_col(
  names(raw),
  c("c_index", "cindex", "c.index", "c-index", "concordance", "concordance_index")
)
split_col <- find_col(names(raw), c("split", "iteration", "repeat", "run", "seed"))
cohort_col <- find_col(names(raw), c("cohort", "dataset", "set"))

# -----------------------------------------------------------------------------
# 4. CHUYEN VE LONG FORMAT
# -----------------------------------------------------------------------------

if (!is.na(model_col) && !is.na(cindex_col)) {
  plot_data <- data.frame(
    split = if (!is.na(split_col)) raw[[split_col]] else seq_len(nrow(raw)),
    cohort = if (!is.na(cohort_col)) as.character(raw[[cohort_col]]) else "TCGA",
    model = canonical_model(raw[[model_col]]),
    c_index = suppressWarnings(as.numeric(raw[[cindex_col]])),
    stringsAsFactors = FALSE
  )
} else {
  # Ho tro wide format: moi model la mot cot.
  mapped_names <- canonical_model(names(raw))
  model_cols <- which(mapped_names %in% MODEL_LEVELS)

  if (length(model_cols) == 0L) {
    stop(
      "Khong nhan dien duoc cot model/C-index trong file:\n", INPUT_FILE,
      "\n\nColumns:\n- ", paste(names(raw), collapse = "\n- ")
    )
  }

  split_values <- if (!is.na(split_col)) raw[[split_col]] else seq_len(nrow(raw))

  plot_data <- rbindlist(lapply(model_cols, function(j) {
    data.frame(
      split = split_values,
      cohort = "TCGA",
      model = mapped_names[j],
      c_index = suppressWarnings(as.numeric(raw[[j]])),
      stringsAsFactors = FALSE
    )
  }), fill = TRUE)
}

plot_data <- as.data.frame(plot_data)

cat("\nCanonical model names detected:\n")
print(unique(plot_data$model))

plot_data <- plot_data[
  is.finite(plot_data$c_index) & plot_data$model %in% MODEL_LEVELS,
  , drop = FALSE
]

if (nrow(plot_data) == 0L) {
  stop("File da doc nhung khong co C-index hop le cho cac model Figure 2A.")
}

# Neu file co ca training/test, uu tien test/validation neu co cot cohort/set.
cohort_key <- key(plot_data$cohort)
if (any(grepl("test|testing|validation|holdout", cohort_key))) {
  keep <- grepl("test|testing|validation|holdout", cohort_key)
  plot_data <- plot_data[keep, , drop = FALSE]
}

# Loai dong trung lap do output co the duoc ghep lai nhieu lan.
plot_data <- unique(plot_data[, c("split", "model", "c_index")])

available_models <- MODEL_LEVELS[MODEL_LEVELS %in% unique(plot_data$model)]
missing_models <- setdiff(MODEL_LEVELS, available_models)

cat("\nModels available (", length(available_models), "/9):\n- ",
    paste(available_models, collapse = "\n- "), "\n", sep = "")

if (length(missing_models) > 0L) {
  cat("\nModels NOT trained in this result file:\n- ",
      paste(missing_models, collapse = "\n- "), "\n", sep = "")
}

# Ghi audit, khong dung stop neu chi co 3 models.
audit <- data.frame(
  model = MODEL_LEVELS,
  status = ifelse(MODEL_LEVELS %in% available_models, "available", "missing"),
  n_values = vapply(MODEL_LEVELS, function(m) sum(plot_data$model == m), integer(1)),
  stringsAsFactors = FALSE
)

fwrite(audit, file.path(RESULT_DIR, "08_Figure2A_model_audit.csv"))
fwrite(plot_data, file.path(RESULT_DIR, "08_Figure2A_plot_data.csv"))

# -----------------------------------------------------------------------------
# 5. SUMMARY
# -----------------------------------------------------------------------------

summary_dt <- as.data.table(plot_data)[,
  .(
    n = .N,
    median_Cindex = median(c_index, na.rm = TRUE),
    mean_Cindex = mean(c_index, na.rm = TRUE),
    Q1 = quantile(c_index, 0.25, na.rm = TRUE),
    Q3 = quantile(c_index, 0.75, na.rm = TRUE),
    min = min(c_index, na.rm = TRUE),
    max = max(c_index, na.rm = TRUE)
  ),
  by = model
]

summary_dt[, order_id := match(model, MODEL_LEVELS)]
setorder(summary_dt, order_id)
summary_dt[, order_id := NULL]
fwrite(summary_dt, file.path(RESULT_DIR, "08_Figure2A_summary.csv"))
print(summary_dt)

# -----------------------------------------------------------------------------
# 6. VE CAC MODEL THUC SU CO KET QUA
# Khong tao gia 6 model con thieu.
# -----------------------------------------------------------------------------

plot_data$model <- factor(plot_data$model, levels = available_models)

model_colors <- c(
  "Clin" = "#BDBDBD",
  "M2EFM Meth+Exp+Clin" = "#E69F00",
  "M2EFM Exp+Clin" = "#56B4E9",
  "Cox Meth+Exp+Clin" = "#009E73",
  "rorS+Clin" = "#F0E442",
  "M2EFM Meth+Exp" = "#D89000",
  "M2EFM Exp" = "#63B8E5",
  "Cox Meth+Exp" = "#159A83",
  "rorS" = "#E8DA3D"
)

subtitle_text <- if (length(available_models) == 9L) {
  "100 random 70/30 splits - full 9-model comparison"
} else {
  paste0(
    "100 random 70/30 splits - ", length(available_models),
    " model(s) available in 03_internal_full_monte_carlo_metrics.csv"
  )
}

p <- ggplot(plot_data, aes(x = model, y = c_index, fill = model)) +
  geom_hline(yintercept = 0.70, linetype = "dashed", linewidth = 0.55, colour = "grey35") +
  geom_boxplot(
    width = 0.70,
    linewidth = 0.55,
    outlier.shape = 16,
    outlier.size = 1.8,
    colour = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 19,
    size = 2.4,
    colour = "black"
  ) +
  scale_fill_manual(values = model_colors[available_models], drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  coord_cartesian(ylim = c(0.43, 0.97)) +
  labs(
    title = "A) TCGA internal validation",
    subtitle = subtitle_text,
    x = NULL,
    y = "C-index",
    fill = "Model"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p)

png_file <- file.path(RESULT_DIR, "08_Figure2A_internal_Cindex.png")
pdf_file <- file.path(RESULT_DIR, "08_Figure2A_internal_Cindex.pdf")

ggsave(png_file, p, width = 12.5, height = 7.2, units = "in", dpi = 300, bg = "white")
ggsave(pdf_file, p, width = 12.5, height = 7.2, units = "in", bg = "white")

cat("\nDONE\n")
cat("Input :", INPUT_FILE, "\n")
cat("PNG   :", png_file, "\n")
cat("PDF   :", pdf_file, "\n")
cat("Audit :", file.path(RESULT_DIR, "08_Figure2A_model_audit.csv"), "\n")

if (length(missing_models) > 0L) {
  cat(
    "\nIMPORTANT: File buoc 03 hien chi co ", length(available_models),
    " model(s). De co 9 cot that, can train them trong 03_train_internal_full.R; ",
    "khong the tao 6 C-index con thieu chi bang code ve.\n",
    sep = ""
  )
}
