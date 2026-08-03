# ============================================================
# 09_plot_Figure2A_9models.R
#
# Draw Figure 2A from:
# results/03_internal_9models_monte_carlo_metrics.csv
# ============================================================

source("config.R")

required_packages <- c("data.table", "ggplot2")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

INPUT_FILE <- file.path(
  RESULTS_DIR,
  "03_internal_9models_monte_carlo_metrics.csv"
)

OUTPUT_PNG <- file.path(
  RESULTS_DIR,
  "09_Figure2A_internal_Cindex_9models.png"
)

OUTPUT_PDF <- file.path(
  RESULTS_DIR,
  "09_Figure2A_internal_Cindex_9models.pdf"
)

OUTPUT_SUMMARY <- file.path(
  RESULTS_DIR,
  "09_Figure2A_internal_Cindex_9models_summary.csv"
)

if (!file.exists(INPUT_FILE)) {
  stop(
    "Cannot find: ", INPUT_FILE,
    "\nRun source(\"03_train_internal_9models.R\") first."
  )
}

x <- data.table::fread(INPUT_FILE)

required_columns <- c(
  "split", "cohort", "model", "c_index"
)

absent <- setdiff(
  required_columns,
  colnames(x)
)

if (length(absent) > 0L) {
  stop(
    "Input file is missing columns: ",
    paste(absent, collapse = ", ")
  )
}

model_levels <- c(
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

x <- x[
  cohort == "TCGA_test" &
    model %in% model_levels &
    is.finite(c_index)
]

counts <- x[, .N, by = model]

missing_models <- setdiff(
  model_levels,
  counts$model
)

wrong_counts <- counts[N != N_MONTE_CARLO_SPLITS]

if (
  length(missing_models) > 0L ||
    nrow(wrong_counts) > 0L
) {
  stop(
    "The 9-model metrics are incomplete.",
    if (length(missing_models) > 0L) {
      paste0(
        "\nMissing models:\n- ",
        paste(missing_models, collapse = "\n- ")
      )
    } else {
      ""
    },
    if (nrow(wrong_counts) > 0L) {
      paste0(
        "\nModels without ",
        N_MONTE_CARLO_SPLITS,
        " rows:\n",
        paste(
          wrong_counts$model,
          wrong_counts$N,
          sep = ": ",
          collapse = "\n"
        )
      )
    } else {
      ""
    }
  )
}

x[, model := factor(
  model,
  levels = model_levels,
  ordered = TRUE
)]

summary_table <- x[
  ,
  .(
    n = .N,
    median_Cindex = stats::median(c_index),
    mean_Cindex = mean(c_index),
    Q1 = stats::quantile(c_index, 0.25),
    Q3 = stats::quantile(c_index, 0.75),
    min = min(c_index),
    max = max(c_index)
  ),
  by = model
]

data.table::fwrite(
  summary_table,
  OUTPUT_SUMMARY
)

model_colors <- c(
  "Clin" = "#BDBDBD",
  "M2EFM Meth+Exp+Clin" = "#E69F00",
  "M2EFM Exp+Clin" = "#56B4E9",
  "Cox Meth+Exp+Clin" = "#009E73",
  "rorS+Clin" = "#F0E442",
  "M2EFM Meth+Exp" = "#E69F00",
  "M2EFM Exp" = "#56B4E9",
  "Cox Meth+Exp" = "#009E73",
  "rorS" = "#F0E442"
)

p <- ggplot2::ggplot(
  x,
  ggplot2::aes(
    x = model,
    y = c_index,
    fill = model
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.70,
    colour = "red",
    linewidth = 0.55
  ) +
  ggplot2::geom_boxplot(
    width = 0.72,
    notch = TRUE,
    colour = "black",
    linewidth = 0.45,
    outlier.shape = 16,
    outlier.size = 1.6,
    outlier.alpha = 0.85
  ) +
  ggplot2::scale_fill_manual(
    values = model_colors,
    limits = model_levels,
    drop = FALSE
  ) +
  ggplot2::scale_x_discrete(
    limits = model_levels,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0.50, 0.90, by = 0.10),
    minor_breaks = seq(0.45, 0.90, by = 0.05)
  ) +
  ggplot2::coord_cartesian(
    ylim = c(0.43, 0.92)
  ) +
  ggplot2::labs(
    tag = "A",
    x = "Model",
    y = "C-index",
    fill = "Model"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 14
    ),
    axis.title = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 48,
      hjust = 1,
      vjust = 1
    ),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(
      face = "bold"
    )
  )

print(p)

ggplot2::ggsave(
  OUTPUT_PNG,
  p,
  width = 13,
  height = 6.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  OUTPUT_PDF,
  p,
  width = 13,
  height = 6.8,
  units = "in",
  bg = "white"
)

message("DONE")
message("Input  : ", INPUT_FILE)
message("PNG    : ", OUTPUT_PNG)
message("PDF    : ", OUTPUT_PDF)
message("Summary: ", OUTPUT_SUMMARY)

print(summary_table)
