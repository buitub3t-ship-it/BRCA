#!/usr/bin/env Rscript

# ============================================================
# STEP 13C — CORRECT FIGURE 3 TIME UNITS AND REDRAW KM PLOTS
#
# Correct source units:
#   TCGA      : days
#   Terunuma  : months
#   Kao       : years
#
# Conversions:
#   Terunuma × (365.25 / 12)
#   Kao       × 365.25
#
# This script does NOT refit the M2EFM model and does NOT change
# risk scores or Low/Medium/High assignments.
#
# Input:
#   results/13_figure3_risk_scores.csv
#
# Outputs:
#   results/13C_figure3_KM_time_units_corrected.png
#   results/13C_figure3_KM_time_units_corrected.pdf
#   results/13C_figure3_time_unit_audit.csv
#   results/13C_figure3_risk_scores_time_corrected.csv
#   results/13C_figure3_logrank_tests.csv
#   results/13C_figure3_hazard_ratios.csv
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1,
  width = 180
)

root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

outdir <- file.path(
  root,
  "results"
)

input_file <- file.path(
  outdir,
  "13_figure3_risk_scores.csv"
)

if (!requireNamespace(
  "survival",
  quietly = TRUE
)) {
  stop(
    "Package 'survival' is required."
  )
}

suppressPackageStartupMessages(
  library(survival)
)

if (!file.exists(
  input_file
)) {
  stop(
    "Missing input: ",
    input_file
  )
}


# ------------------------------------------------------------
# 1. Fixed source-specific time units
# ------------------------------------------------------------

cohort_order <- c(
  "TCGA",
  "Terunuma",
  "Kao"
)

risk_levels <- c(
  "Low",
  "Medium",
  "High"
)

panel_letters <- c(
  TCGA = "A",
  Terunuma = "B",
  Kao = "C"
)

risk_colors <- c(
  Low = "#00008B",
  Medium = "#FF9800",
  High = "#FF0000"
)

source_units <- c(
  TCGA = "days",
  Terunuma = "months",
  Kao = "years"
)

multiplier_to_days <- c(
  TCGA = 1,
  Terunuma = 365.25 / 12,
  Kao = 365.25
)

plot_limit_days <- 3652.5

x_breaks <- c(
  0,
  1000,
  2000,
  3000
)


# ------------------------------------------------------------
# 2. Read Step 13 risk scores and validate
# ------------------------------------------------------------

x <- read.csv(
  input_file,
  check.names = FALSE
)

required_columns <- c(
  "sample_id",
  "cohort",
  "time_days",
  "event",
  "risk_group"
)

missing_columns <- setdiff(
  required_columns,
  colnames(x)
)

if (length(
  missing_columns
) > 0L) {
  stop(
    "Missing columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

x$cohort <- as.character(
  x$cohort
)

x$time_original <- suppressWarnings(
  as.numeric(
    x$time_days
  )
)

x$event <- suppressWarnings(
  as.numeric(
    x$event
  )
)

x$risk_group <- factor(
  as.character(
    x$risk_group
  ),
  levels = risk_levels
)

if (
  anyNA(
    x$time_original
  ) ||
    anyNA(
      x$event
    ) ||
    anyNA(
      x$risk_group
    ) ||
    any(
      !is.finite(
        x$time_original
      )
    ) ||
    any(
      x$time_original < 0
    ) ||
    any(
      !x$event %in% c(
        0,
        1
      )
    )
) {
  stop(
    "Invalid survival, event, or risk-group values."
  )
}

if (!setequal(
  unique(
    x$cohort
  ),
  cohort_order
)) {
  stop(
    "Unexpected cohort names: ",
    paste(
      sort(
        unique(
          x$cohort
        )
      ),
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 3. Convert source times to days
# ------------------------------------------------------------

audit_rows <- list()

for (cohort in cohort_order) {

  index <- x$cohort ==
    cohort

  multiplier <- multiplier_to_days[
    [cohort]
  ]

  x$time_days[index] <-
    x$time_original[index] *
    multiplier

  audit_rows[[cohort]] <- data.frame(
    cohort = cohort,
    source_unit = source_units[
      [cohort]
    ],
    multiplier_to_days = multiplier,
    samples = sum(
      index
    ),
    original_min = min(
      x$time_original[index]
    ),
    original_median = median(
      x$time_original[index]
    ),
    original_max = max(
      x$time_original[index]
    ),
    corrected_min_days = min(
      x$time_days[index]
    ),
    corrected_median_days = median(
      x$time_days[index]
    ),
    corrected_max_days = max(
      x$time_days[index]
    ),
    stringsAsFactors = FALSE
  )
}

time_audit <- do.call(
  rbind,
  audit_rows
)

rownames(
  time_audit
) <- NULL

if (
  anyNA(
    x$time_days
  ) ||
    any(
      !is.finite(
        x$time_days
      )
    )
) {
  stop(
    "Time conversion produced invalid values."
  )
}

write.csv(
  time_audit,
  file.path(
    outdir,
    "13C_figure3_time_unit_audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  x,
  file.path(
    outdir,
    "13C_figure3_risk_scores_time_corrected.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. KM estimates, log-rank tests, and risk-group HRs
# ------------------------------------------------------------

km_fits <- setNames(
  vector(
    "list",
    length(
      cohort_order
    )
  ),
  cohort_order
)

logrank_rows <- list()
hazard_ratio_rows <- list()

for (cohort in cohort_order) {

  dat <- x[
    x$cohort ==
      cohort,
    ,
    drop = FALSE
  ]

  dat$risk_group <- factor(
    dat$risk_group,
    levels = risk_levels
  )

  km_fits[[cohort]] <- survival::survfit(
    survival::Surv(
      time_days,
      event
    ) ~ risk_group,
    data = dat
  )

  logrank_fit <- survival::survdiff(
    survival::Surv(
      time_days,
      event
    ) ~ risk_group,
    data = dat,
    rho = 0
  )

  logrank_df <- length(
    logrank_fit$n
  ) - 1L

  logrank_rows[[cohort]] <- data.frame(
    cohort = cohort,
    chi_square = unname(
      logrank_fit$chisq
    ),
    degrees_freedom = logrank_df,
    p_value = stats::pchisq(
      logrank_fit$chisq,
      df = logrank_df,
      lower.tail = FALSE
    ),
    stringsAsFactors = FALSE
  )

  cox_fit <- survival::coxph(
    survival::Surv(
      time_days,
      event
    ) ~ risk_group,
    data = dat,
    ties = "breslow"
  )

  cox_summary <- summary(
    cox_fit
  )

  terms <- c(
    "risk_groupMedium",
    "risk_groupHigh"
  )

  comparisons <- c(
    "Medium vs Low",
    "High vs Low"
  )

  cohort_hr_rows <- lapply(
    seq_along(
      terms
    ),
    function(i) {

      term <- terms[
        [i]
      ]

      data.frame(
        cohort = cohort,
        comparison = comparisons[
          [i]
        ],
        hazard_ratio = cox_summary$conf.int[
          term,
          "exp(coef)"
        ],
        lower_95 = cox_summary$conf.int[
          term,
          "lower .95"
        ],
        upper_95 = cox_summary$conf.int[
          term,
          "upper .95"
        ],
        p_value = cox_summary$coefficients[
          term,
          "Pr(>|z|)"
        ],
        stringsAsFactors = FALSE
      )
    }
  )

  hazard_ratio_rows[[cohort]] <- do.call(
    rbind,
    cohort_hr_rows
  )
}

logrank_table <- do.call(
  rbind,
  logrank_rows
)

rownames(
  logrank_table
) <- NULL

hazard_ratio_table <- do.call(
  rbind,
  hazard_ratio_rows
)

rownames(
  hazard_ratio_table
) <- NULL

write.csv(
  logrank_table,
  file.path(
    outdir,
    "13C_figure3_logrank_tests.csv"
  ),
  row.names = FALSE
)

write.csv(
  hazard_ratio_table,
  file.path(
    outdir,
    "13C_figure3_hazard_ratios.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Author-style KM figure with full legend
# ------------------------------------------------------------

draw_figure <- function(
    output_file,
    device = c(
      "png",
      "pdf"
    )
) {

  device <- match.arg(
    device
  )

  if (device ==
    "png") {
    grDevices::png(
      filename = output_file,
      width = 4200,
      height = 1450,
      res = 320,
      bg = "white"
    )
  } else {
    grDevices::pdf(
      file = output_file,
      width = 13.1,
      height = 4.55,
      onefile = TRUE
    )
  }

  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )

  old_par <- graphics::par(
    no.readonly = TRUE
  )

  on.exit(
    graphics::par(
      old_par
    ),
    add = TRUE
  )

  graphics::layout(
    matrix(
      1:4,
      nrow = 1
    ),
    widths = c(
      1,
      1,
      1,
      0.95
    )
  )

  for (cohort in cohort_order) {

    graphics::par(
      mar = c(
        5.0,
        4.5,
        2.2,
        1.0
      ),
      mgp = c(
        2.65,
        0.78,
        0
      ),
      tcl = -0.28,
      las = 1,
      bty = "l",
      xpd = FALSE
    )

    plot(
      km_fits[[cohort]],
      col = unname(
        risk_colors
      ),
      lwd = 1.65,
      lty = 1,
      mark.time = FALSE,
      conf.int = FALSE,
      xlim = c(
        0,
        plot_limit_days
      ),
      ylim = c(
        0.25,
        1.02
      ),
      xlab = "Days",
      ylab = "Survival",
      xaxt = "n",
      yaxt = "n",
      cex.lab = 1.0
    )

    graphics::axis(
      side = 1,
      at = x_breaks,
      labels = as.character(
        x_breaks
      ),
      las = 2,
      cex.axis = 0.86
    )

    graphics::axis(
      side = 2,
      at = seq(
        0.4,
        1.0,
        by = 0.2
      ),
      labels = format(
        seq(
          0.4,
          1.0,
          by = 0.2
        ),
        nsmall = 1
      ),
      las = 1,
      cex.axis = 0.86
    )

    graphics::mtext(
      panel_letters[
        [cohort]
      ],
      side = 3,
      line = 0.45,
      adj = -0.08,
      font = 2,
      cex = 1.2
    )
  }

  graphics::par(
    mar = c(
      0,
      0.4,
      0,
      0.5
    ),
    xpd = NA
  )

  graphics::plot.new()

  graphics::legend(
    x = "center",
    legend = risk_levels,
    title = "Risk Group",
    col = unname(
      risk_colors
    ),
    lty = 1,
    lwd = 1.65,
    bty = "n",
    seg.len = 2.0,
    cex = 1.02,
    x.intersp = 0.8,
    y.intersp = 1.1
  )

  invisible(
    NULL
  )
}

png_file <- file.path(
  outdir,
  "13C_figure3_KM_time_units_corrected.png"
)

pdf_file <- file.path(
  outdir,
  "13C_figure3_KM_time_units_corrected.pdf"
)

draw_figure(
  png_file,
  "png"
)

draw_figure(
  pdf_file,
  "pdf"
)


# ------------------------------------------------------------
# 6. Console report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 13C COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nTIME-UNIT AUDIT\n")
print(
  time_audit,
  row.names = FALSE
)

cat("\nLOG-RANK TESTS\n")
print(
  logrank_table,
  row.names = FALSE
)

cat("\nHAZARD RATIOS\n")
print(
  hazard_ratio_table,
  row.names = FALSE
)

cat("\nCorrected figure saved:\n")
cat(
  png_file,
  "\n",
  pdf_file,
  "\n",
  sep = ""
)

cat("\nFinished successfully.\n")
