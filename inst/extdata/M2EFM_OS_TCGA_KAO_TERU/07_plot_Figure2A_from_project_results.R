# ============================================================================
# 07_plot_Figure2A_from_project_results.R
#
# Muc dich:
#   - KHONG yeu cau file gia dinh results/03_internal_Cindex_9models.csv.
#   - Tu dong tim cac file CSV trong results/ va csv_output/.
#   - Chuan hoa ten mo hinh va ve Figure 2A khi da co DU 9 mo hinh.
#   - Khong tu tao/gia lap C-index cho cac mo hinh chua duoc train.
#
# Chay trong RStudio project M2EFM_OS_TCGA_KAO_TERU:
#   source("07_plot_Figure2A_from_project_results.R")
# ============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Thieu package: ", paste(missing_packages, collapse = ", "),
    "\nCai bang install.packages(c(\"data.table\", \"ggplot2\"))."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# 1. PATH INPUT / OUTPUT
# ---------------------------------------------------------------------------

# Thu muc project dang mo trong anh cua ban:
PROJECT_DIR <- normalizePath(
  "D:/M2EFM_BRCA_Replication/M2EFM-master/inst/extdata/M2EFM_OS_TCGA_KAO_TERU",
  winslash = "/",
  mustWork = TRUE
)

# Script se tim CSV ket qua trong ca hai thu muc nay.
# - results/: C-index va cac output tu 03/04/06
# - csv_output/: cac CSV trong project; file raw khong co c_index se tu dong bo qua
INPUT_DIRS <- c(
  file.path(PROJECT_DIR, "results"),
  file.path(PROJECT_DIR, "csv_output")
)
INPUT_DIRS <- INPUT_DIRS[dir.exists(INPUT_DIRS)]

if (length(INPUT_DIRS) == 0L) {
  stop(
    "Khong tim thay thu muc input. Da kiem tra:\n",
    paste(c(
      file.path(PROJECT_DIR, "results"),
      file.path(PROJECT_DIR, "csv_output")
    ), collapse = "\n")
  )
}

# Tat ca output cua script duoc luu tai day.
OUTPUT_DIR <- file.path(PROJECT_DIR, "results", "Figure2A_9models")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Project directory :", PROJECT_DIR, "\n")
cat("Input directories:", paste(INPUT_DIRS, collapse = " ; "), "\n")
cat("Output directory :", OUTPUT_DIR, "\n\n")

# ---------------------------------------------------------------------------
# 2. DANH SACH 9 MO HINH FIGURE 2A
# ---------------------------------------------------------------------------

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

normalize_key <- function(x) {
  tolower(gsub("[^a-z0-9]", "", as.character(x)))
}

# Chuan hoa ten mo hinh. Hai ten generic cua project cu duoc hieu theo context:
# - file internal/full/meta: M2EFM_ExpClin = Meth+Exp+Clin;
#                            M2EFM_Molecular = Meth+Exp.
# - file expression/external: M2EFM_ExpClin = Exp+Clin;
#                             M2EFM_Molecular = Exp.
canonical_model <- function(model, source_file = "") {
  key <- normalize_key(model)
  source_key <- normalize_key(source_file)

  integrated_context <- grepl(
    "figure2a|9models|9model|internalfull|metadimensional|meta|methyl|methexp|03traininternal",
    source_key
  ) && !grepl("externalexpression|04trainexternal", source_key)

  out <- rep(NA_character_, length(key))

  out[key %in% c("clin", "clinical", "clinicalonly", "covariates")] <- "Clin"

  out[key %in% c(
    "m2efmmethexpclin", "m2efmmethylationexpressionclinical",
    "m2efmintegratedclin", "m2efmintegratedclinical"
  )] <- "M2EFM Meth+Exp+Clin"

  out[key %in% c(
    "m2efmexpclin", "m2efmexpressionclin", "m2efmexpressionclinical"
  )] <- "M2EFM Exp+Clin"

  out[key %in% c(
    "coxmethexpclin", "coxridgemethexpclin", "coxintegratedclin",
    "coxridgemethylationexpressionclinical"
  )] <- "Cox Meth+Exp+Clin"

  out[key %in% c("rorsclin", "pam50rorsclin", "pam50riskclin")] <- "rorS+Clin"

  out[key %in% c(
    "m2efmmethexp", "m2efmmethylationexpression", "m2efmintegrated"
  )] <- "M2EFM Meth+Exp"

  out[key %in% c("m2efmexp", "m2efmexpression")] <- "M2EFM Exp"

  out[key %in% c(
    "coxmethexp", "coxridgemethexp", "coxintegrated",
    "coxridgemethylationexpression"
  )] <- "Cox Meth+Exp"

  out[key %in% c("rors", "pam50rors", "pam50risk")] <- "rorS"

  # Ten generic cua cac script da chay truoc day.
  generic_expclin <- key == "m2efmexpclin"
  generic_mol     <- key == "m2efmmolecular"

  if (integrated_context) {
    out[generic_expclin] <- "M2EFM Meth+Exp+Clin"
    out[generic_mol]     <- "M2EFM Meth+Exp"
  } else {
    out[generic_expclin] <- "M2EFM Exp+Clin"
    out[generic_mol]     <- "M2EFM Exp"
  }

  out
}

# ---------------------------------------------------------------------------
# 3. DOC MOT FILE CSV KET QUA
# ---------------------------------------------------------------------------

find_column <- function(nms, candidates) {
  nkey <- normalize_key(nms)
  ckey <- normalize_key(candidates)
  hit <- match(ckey, nkey)
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0L) return(NA_character_)
  nms[hit[1L]]
}

read_result_csv <- function(path) {
  z <- tryCatch(
    data.table::fread(
      path,
      data.table = FALSE,
      na.strings = c("NA", "N/A", "--", "")
    ),
    error = function(e) NULL
  )

  if (is.null(z) || nrow(z) == 0L || ncol(z) == 0L) return(NULL)

  model_col <- find_column(names(z), c("model", "method", "model_name"))
  cindex_col <- find_column(
    names(z),
    c("c_index", "cindex", "c-index", "concordance", "concordance_index")
  )
  split_col <- find_column(names(z), c("split", "iteration", "repeat", "run", "seed"))
  cohort_col <- find_column(names(z), c("cohort", "dataset", "set", "sample_set"))

  # Long format: split/cohort/model/c_index
  if (!is.na(model_col) && !is.na(cindex_col)) {
    out <- data.frame(
      source_file = path,
      split = if (!is.na(split_col)) z[[split_col]] else seq_len(nrow(z)),
      cohort = if (!is.na(cohort_col)) as.character(z[[cohort_col]]) else NA_character_,
      raw_model = as.character(z[[model_col]]),
      c_index = suppressWarnings(as.numeric(z[[cindex_col]])),
      stringsAsFactors = FALSE
    )

    out$model <- canonical_model(out$raw_model, path)
    return(out)
  }

  # Wide format: moi cot la mot mo hinh.
  mapped <- canonical_model(names(z), path)
  model_columns <- which(!is.na(mapped))

  if (length(model_columns) == 0L) return(NULL)

  id_split <- if (!is.na(split_col)) z[[split_col]] else seq_len(nrow(z))
  id_cohort <- if (!is.na(cohort_col)) as.character(z[[cohort_col]]) else NA_character_

  out_list <- lapply(model_columns, function(j) {
    data.frame(
      source_file = path,
      split = id_split,
      cohort = id_cohort,
      raw_model = names(z)[j],
      model = mapped[j],
      c_index = suppressWarnings(as.numeric(z[[j]])),
      stringsAsFactors = FALSE
    )
  })

  data.table::rbindlist(out_list, fill = TRUE)
}

# ---------------------------------------------------------------------------
# 4. QUET TOAN BO RESULTS/
# ---------------------------------------------------------------------------

csv_files <- unique(unlist(lapply(INPUT_DIRS, function(input_dir) {
  list.files(
    input_dir,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
})))

# Khong doc lai output cua chinh script nay khi chay lan sau.
csv_files <- csv_files[
  !grepl("Figure2A_9models", csv_files, fixed = TRUE)
]

if (length(csv_files) == 0L) {
  stop("Khong tim thay file CSV nao trong: ", paste(INPUT_DIRS, collapse = "; "))
}

cat("CSV files found:", length(csv_files), "\n")

parts <- lapply(csv_files, read_result_csv)
parts <- parts[!vapply(parts, is.null, logical(1))]

if (length(parts) == 0L) {
  stop(
    "Da tim thay CSV, nhung khong file nao co cot model + c_index ",
    "hoac cac cot wide mang ten mo hinh."
  )
}

all_rows <- data.table::rbindlist(parts, fill = TRUE)
all_rows <- as.data.frame(all_rows)

all_rows <- all_rows[
  !is.na(all_rows$model) &
    is.finite(all_rows$c_index) &
    all_rows$c_index >= 0 &
    all_rows$c_index <= 1,
  ,
  drop = FALSE
]

# Figure 2A la TCGA internal test. Loai cac cohort external ro rang.
cohort_key <- normalize_key(all_rows$cohort)
external_cohort <- cohort_key %in% c(
  "kao", "terunuma", "hatzis1", "hatzis2", "validation", "external"
)
all_rows <- all_rows[is.na(external_cohort) | !external_cohort, , drop = FALSE]

if (nrow(all_rows) == 0L) {
  stop("Khong con dong C-index TCGA/internal hop le sau khi loc.")
}

# Loai trung lap trong cung mot file.
all_rows <- unique(all_rows[, c(
  "source_file", "split", "cohort", "raw_model", "model", "c_index"
)])

# Audit tat ca dong ma script tim thay.
write.csv(
  all_rows,
  file.path(OUTPUT_DIR, "07_discovered_model_rows.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# 5. CHON NGUON TOT NHAT CHO TUNG MO HINH
# ---------------------------------------------------------------------------

source_stats <- aggregate(
  c_index ~ model + source_file,
  data = all_rows,
  FUN = length
)
names(source_stats)[names(source_stats) == "c_index"] <- "n"

source_key <- normalize_key(source_stats$source_file)
source_stats$priority <- 0L
source_stats$priority <- source_stats$priority +
  ifelse(grepl("figure2a|9models|9model", source_key), 1000L, 0L)
source_stats$priority <- source_stats$priority +
  ifelse(grepl("internalfull|03traininternal|internal", source_key), 100L, 0L)
source_stats$priority <- source_stats$priority -
  ifelse(grepl("external|04trainexternal|terunuma|kao", source_key), 500L, 0L)
source_stats$score <- source_stats$priority + pmin(source_stats$n, 200L)

source_stats <- source_stats[order(
  source_stats$model,
  -source_stats$score,
  -source_stats$n,
  source_stats$source_file
), ]

best_sources <- source_stats[!duplicated(source_stats$model), ]

plot_rows <- merge(
  all_rows,
  best_sources[, c("model", "source_file")],
  by = c("model", "source_file")
)

# Mot so file summary chi co 1 median/model. Khong duoc dung de tao boxplot.
counts <- aggregate(c_index ~ model, data = plot_rows, FUN = length)
names(counts)[2] <- "n"

write.csv(
  best_sources,
  file.path(OUTPUT_DIR, "07_selected_source_per_model.csv"),
  row.names = FALSE
)

write.csv(
  counts,
  file.path(OUTPUT_DIR, "07_model_counts.csv"),
  row.names = FALSE
)

cat("\nModels discovered:\n")
print(counts, row.names = FALSE)

available <- unique(plot_rows$model)
missing_models <- setdiff(MODEL_LEVELS, available)
insufficient <- counts$model[counts$n < 2L]

if (length(missing_models) > 0L || length(insufficient) > 0L) {
  message("\n============================================================")
  message("CHUA THE VE DU 9 BOXPLOT")

  if (length(missing_models) > 0L) {
    message("Cac mo hinh chua co ket qua per-split:")
    message("- ", paste(missing_models, collapse = "\n- "))
  }

  if (length(insufficient) > 0L) {
    message("Cac mo hinh chi co summary/qua it dong:")
    message("- ", paste(insufficient, collapse = "\n- "))
  }

  message("\nAudit da luu tai:")
  message(file.path(OUTPUT_DIR, "07_discovered_model_rows.csv"))
  message(file.path(OUTPUT_DIR, "07_selected_source_per_model.csv"))
  message("============================================================")

  stop(
    "Project hien tai chua co du ket qua C-index theo tung split cho 9 mo hinh. ",
    "Day la thieu ket qua training, khong phai loi duong dan cua plotting script."
  )
}

# ---------------------------------------------------------------------------
# 6. SUMMARY + PLOT
# ---------------------------------------------------------------------------

plot_rows$model <- factor(
  plot_rows$model,
  levels = MODEL_LEVELS,
  ordered = TRUE
)

summary_table <- do.call(
  rbind,
  lapply(split(plot_rows$c_index, plot_rows$model), function(x) {
    data.frame(
      n = length(x),
      median = median(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      q025 = unname(quantile(x, 0.025, na.rm = TRUE)),
      q25 = unname(quantile(x, 0.25, na.rm = TRUE)),
      q75 = unname(quantile(x, 0.75, na.rm = TRUE)),
      q975 = unname(quantile(x, 0.975, na.rm = TRUE))
    )
  })
)
summary_table$model <- rownames(summary_table)
rownames(summary_table) <- NULL
summary_table <- summary_table[, c("model", setdiff(names(summary_table), "model"))]

write.csv(
  plot_rows,
  file.path(OUTPUT_DIR, "07_Figure2A_9models_long.csv"),
  row.names = FALSE
)
write.csv(
  summary_table,
  file.path(OUTPUT_DIR, "07_Figure2A_9models_summary.csv"),
  row.names = FALSE
)

model_colours <- c(
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

p <- ggplot(plot_rows, aes(x = model, y = c_index, fill = model)) +
  geom_hline(yintercept = 0.70, colour = "red", linewidth = 0.55) +
  geom_boxplot(
    width = 0.72,
    linewidth = 0.45,
    colour = "black",
    outlier.shape = 16,
    outlier.size = 1.6,
    outlier.alpha = 0.75
  ) +
  scale_fill_manual(values = model_colours, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = seq(0.4, 0.9, 0.1)) +
  coord_cartesian(ylim = c(0.40, 0.92)) +
  labs(
    tag = "A",
    title = "TCGA internal validation",
    subtitle = paste0("Random 70/30 splits; n per model = ", min(counts$n)),
    x = NULL,
    y = "C-index"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.tag = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(p)

ggsave(
  file.path(OUTPUT_DIR, "07_Figure2A_9models.png"),
  p,
  width = 13,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(OUTPUT_DIR, "07_Figure2A_9models.pdf"),
  p,
  width = 13,
  height = 7,
  units = "in",
  bg = "white"
)

cat("\nCompleted. Outputs:\n")
cat(file.path(OUTPUT_DIR, "07_Figure2A_9models.png"), "\n")
cat(file.path(OUTPUT_DIR, "07_Figure2A_9models.pdf"), "\n")
cat(file.path(OUTPUT_DIR, "07_Figure2A_9models_long.csv"), "\n")
cat(file.path(OUTPUT_DIR, "07_Figure2A_9models_summary.csv"), "\n")
