#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
outdir <- file.path(root, "results")
input_file <- file.path(outdir, "13_figure3_risk_scores.csv")

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package survival is required.")
}
if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

cohort_order <- c("TCGA", "Terunuma", "Kao")
risk_levels <- c("Low", "Medium", "High")
risk_colors <- c(Low = "#00008B", Medium = "#FF9800", High = "#FF0000")
panel_letters <- c(TCGA = "A", Terunuma = "B", Kao = "C")

days_per_month <- 365.25 / 12
plot_limit_days <- 3652.5
x_breaks <- c(0, 1000, 2000, 3000)

x <- read.csv(input_file, check.names = FALSE)

required <- c("sample_id", "cohort", "time_days", "event", "risk_group")
missing <- setdiff(required, colnames(x))
if (length(missing) > 0L) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

x$cohort <- as.character(x$cohort)
x$time_original <- as.numeric(x$time_days)
x$event <- as.numeric(x$event)
x$risk_group <- factor(as.character(x$risk_group), levels = risk_levels)

if (
  anyNA(x$time_original) ||
  anyNA(x$event) ||
  anyNA(x$risk_group) ||
  any(!is.finite(x$time_original)) ||
  any(x$time_original < 0) ||
  any(!x$event %in% c(0, 1))
) {
  stop("Invalid input values.")
}

if (!setequal(unique(x$cohort), cohort_order)) {
  stop("Unexpected cohort names.")
}

detect_unit <- function(values, cohort) {
  maximum <- max(values)
  median_value <- median(values)

  if (cohort == "TCGA") {
    return(list(unit = "days", multiplier = 1, rule = "TCGA retained as days"))
  }

  if (maximum <= 500 && median_value <= 250) {
    return(list(
      unit = "months",
      multiplier = days_per_month,
      rule = "external max<=500 and median<=250"
    ))
  }

  list(
    unit = "days",
    multiplier = 1,
    rule = "already compatible with days"
  )
}

audit_rows <- list()

for (cohort in cohort_order) {
  idx <- x$cohort == cohort
  decision <- detect_unit(x$time_original[idx], cohort)

  x$time_days[idx] <- x$time_original[idx] * decision$multiplier

  audit_rows[[cohort]] <- data.frame(
    cohort = cohort,
    samples = sum(idx),
    original_min = min(x$time_original[idx]),
    original_median = median(x$time_original[idx]),
    original_max = max(x$time_original[idx]),
    detected_unit = decision$unit,
    multiplier_to_days = decision$multiplier,
    corrected_min_days = min(x$time_days[idx]),
    corrected_median_days = median(x$time_days[idx]),
    corrected_max_days = max(x$time_days[idx]),
    detection_rule = decision$rule,
    stringsAsFactors = FALSE
  )
}

time_audit <- do.call(rbind, audit_rows)
rownames(time_audit) <- NULL

write.csv(
  time_audit,
  file.path(outdir, "13B_figure3_time_unit_audit.csv"),
  row.names = FALSE
)

write.csv(
  x,
  file.path(outdir, "13B_figure3_risk_scores_time_corrected.csv"),
  row.names = FALSE
)

km_fits <- setNames(vector("list", length(cohort_order)), cohort_order)

for (cohort in cohort_order) {
  dat <- x[x$cohort == cohort, , drop = FALSE]
  dat$risk_group <- factor(dat$risk_group, levels = risk_levels)

  km_fits[[cohort]] <- survival::survfit(
    survival::Surv(time_days, event) ~ risk_group,
    data = dat
  )
}

draw_figure <- function(file, device = c("png", "pdf")) {
  device <- match.arg(device)

  if (device == "png") {
    png(file, width = 3300, height = 1250, res = 300, bg = "white")
  } else {
    pdf(file, width = 11, height = 4.2, onefile = TRUE)
  }
  on.exit(dev.off(), add = TRUE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  layout(matrix(1:4, nrow = 1), widths = c(1, 1, 1, 0.72))

  for (cohort in cohort_order) {
    par(
      mar = c(4.8, 4.3, 2.0, 0.8),
      mgp = c(2.5, 0.72, 0),
      tcl = -0.28,
      las = 1,
      bty = "l"
    )

    plot(
      km_fits[[cohort]],
      col = unname(risk_colors),
      lwd = 1.55,
      lty = 1,
      mark.time = FALSE,
      conf.int = FALSE,
      xlim = c(0, plot_limit_days),
      ylim = c(0.25, 1.02),
      xlab = "Days",
      ylab = "Survival",
      xaxt = "n",
      yaxt = "n"
    )

    axis(
      1,
      at = x_breaks,
      labels = as.character(x_breaks),
      las = 2,
      cex.axis = 0.88
    )

    axis(
      2,
      at = seq(0.4, 1.0, by = 0.2),
      labels = format(seq(0.4, 1.0, by = 0.2), nsmall = 1),
      las = 1,
      cex.axis = 0.88
    )

    mtext(
      panel_letters[[cohort]],
      side = 3,
      line = 0.45,
      adj = -0.08,
      font = 2,
      cex = 1.18
    )
  }

  par(mar = c(0, 0.5, 0, 0))
  plot.new()

  legend(
    "center",
    legend = risk_levels,
    title = "Risk Group",
    col = unname(risk_colors),
    lty = 1,
    lwd = 1.55,
    bty = "n",
    seg.len = 1.6,
    cex = 1.0
  )
}

png_file <- file.path(outdir, "13B_figure3_KM_axes_fixed.png")
pdf_file <- file.path(outdir, "13B_figure3_KM_axes_fixed.pdf")

draw_figure(png_file, "png")
draw_figure(pdf_file, "pdf")

cat("\n============================================================\n")
cat("STEP 13B COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")
cat("\nTIME-UNIT AUDIT\n")
print(time_audit, row.names = FALSE)
cat("\nCorrected figure saved:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat("\nFinished successfully.\n")
