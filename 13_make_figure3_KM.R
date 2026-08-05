#!/usr/bin/env Rscript

# STEP 13 — Figure 3: final full-TCGA M2EFM Exp+Clin model
# Risk groups follow the paper:
# Low < cohort Q25; Medium Q25–Q75; High > cohort Q75.
# KM curves are displayed through 10 years (3652.5 days).

options(stringsAsFactors = FALSE, warn = 1, width = 180)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
outdir <- file.path(root, "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
input_file <- file.path(outdir, "10B_external_os_combat_rors_ready.rds")

pkgs <- c("glmnet", "survival", "matrixStats", "survcomp")
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Missing packages: ", paste(names(ok)[!ok], collapse = ", "))

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(matrixStats)
  library(survcomp)
})

if (!file.exists(input_file)) stop("Missing input: ", input_file)

cohorts <- c("TCGA", "Terunuma", "Kao")
panel_letters <- c(TCGA = "A", Terunuma = "B", Kao = "C")
risk_levels <- c("Low", "Medium", "High")
risk_colors <- c(Low = "#00008B", Medium = "#FF9800", High = "#FF0000")
stage_levels <- c("Stage I", "Stage II", "Stage III", "Stage IV")

time_col <- "OVERALL.SURVIVAL"
event_col <- "overall.survival.indicator"
age_col <- "age.Dx"
stage_col <- "pathologic_stage"
ten_year_days <- 3652.5
set.seed(1)

validate_expression <- function(x, cohort) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(cohort, ": missing expression names.")
  }
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop(cohort, ": duplicated expression names.")
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(cohort, ": non-finite expression values.")
  }
  x
}

prepare_clinical <- function(x, sample_ids, cohort) {
  req <- c(time_col, event_col, age_col, stage_col)
  miss <- setdiff(req, colnames(x))
  if (length(miss)) stop(cohort, ": missing clinical columns: ", paste(miss, collapse = ", "))
  if (is.null(rownames(x)) || !identical(rownames(x), sample_ids)) {
    stop(cohort, ": expression/clinical samples are not aligned.")
  }

  x[[time_col]] <- as.numeric(x[[time_col]])
  x[[event_col]] <- as.numeric(x[[event_col]])
  x[[age_col]] <- as.numeric(x[[age_col]])
  x[[stage_col]] <- factor(as.character(x[[stage_col]]), levels = stage_levels)

  valid <- complete.cases(x[, req, drop = FALSE]) &
    is.finite(x[[time_col]]) &
    is.finite(x[[event_col]]) &
    is.finite(x[[age_col]]) &
    x[[time_col]] >= 0 &
    x[[event_col]] %in% c(0, 1)

  if (!all(valid)) stop(cohort, ": invalid clinical/survival values.")
  x
}

scale_params <- function(tcga) {
  means <- matrixStats::rowMeans2(tcga)
  ranges <- matrixStats::rowMaxs(tcga) - matrixStats::rowMins(tcga)
  names(means) <- rownames(tcga)
  names(ranges) <- rownames(tcga)
  if (any(!is.finite(ranges) | ranges <= 0)) stop("Invalid TCGA gene ranges.")
  list(means = means, ranges = ranges)
}

apply_scaling <- function(x, params, cohort) {
  if (!identical(rownames(x), names(params$means))) stop(cohort, ": gene order mismatch.")
  z <- sweep(x, 1L, params$means, "-")
  z <- sweep(z, 1L, params$ranges, "/")
  if (anyNA(z) || any(!is.finite(z))) stop(cohort, ": scaling failed.")
  z
}

assign_groups <- function(score) {
  q <- unname(quantile(score, c(0.25, 0.75), type = 7))
  if (any(!is.finite(q)) || q[1] >= q[2]) stop("Invalid risk-score quartiles.")
  group <- ifelse(score < q[1], "Low", ifelse(score > q[2], "High", "Medium"))
  list(group = factor(group, levels = risk_levels), q25 = q[1], q75 = q[2])
}

cindex_row <- function(score, clin, cohort) {
  z <- survcomp::concordance.index(
    x = score,
    surv.time = clin[[time_col]],
    surv.event = clin[[event_col]],
    method = "noether"
  )
  data.frame(
    cohort = cohort,
    c_index = as.numeric(z$c.index),
    lower_95 = as.numeric(z$lower),
    upper_95 = as.numeric(z$upper),
    stringsAsFactors = FALSE
  )
}

hr_rows <- function(dat, cohort) {
  fit <- survival::coxph(
    survival::Surv(time_days, event) ~ risk_group,
    data = dat,
    ties = "breslow"
  )
  s <- summary(fit)
  terms <- c("risk_groupMedium", "risk_groupHigh")
  labels <- c("Medium vs Low", "High vs Low")

  rows <- lapply(seq_along(terms), function(i) {
    term <- terms[i]
    data.frame(
      cohort = cohort,
      comparison = labels[i],
      hazard_ratio = s$conf.int[term, "exp(coef)"],
      lower_95 = s$conf.int[term, "lower .95"],
      upper_95 = s$conf.int[term, "upper .95"],
      p_value = s$coefficients[term, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  })
  list(fit = fit, table = do.call(rbind, rows))
}

logrank_row <- function(dat, cohort) {
  z <- survival::survdiff(
    survival::Surv(time_days, event) ~ risk_group,
    data = dat
  )
  df <- length(z$n) - 1L
  data.frame(
    cohort = cohort,
    chi_square = unname(z$chisq),
    degrees_freedom = df,
    p_value = pchisq(z$chisq, df = df, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

cat("Loading Step 10B object...\n")
ready <- readRDS(input_file)

missing_cohorts <- setdiff(cohorts, names(ready))
if (length(missing_cohorts)) {
  stop("Input lacks cohorts: ", paste(missing_cohorts, collapse = ", "))
}

expression <- setNames(lapply(cohorts, function(cohort) {
  validate_expression(ready[[cohort]]$expression, cohort)
}), cohorts)

clinical <- setNames(lapply(cohorts, function(cohort) {
  prepare_clinical(
    ready[[cohort]]$clinical,
    colnames(expression[[cohort]]),
    cohort
  )
}), cohorts)

common_genes <- rownames(expression$TCGA)
for (cohort in cohorts) {
  if (!identical(rownames(expression[[cohort]]), common_genes)) {
    stop(cohort, ": common-gene order mismatch.")
  }
}

signature <- unique(as.character(ready$external_expression_signature))
if (length(signature) != 113L) {
  stop("Expected 113 common M2EFM genes; found ", length(signature), ".")
}
missing_signature <- setdiff(signature, common_genes)
if (length(missing_signature)) {
  stop("Signature genes missing: ", paste(missing_signature, collapse = ", "))
}

cat("Scaling expression with full-TCGA gene means and ranges...\n")
params <- scale_params(expression$TCGA)
scaled <- setNames(lapply(cohorts, function(cohort) {
  apply_scaling(expression[[cohort]], params, cohort)
}), cohorts)

m2efm_x <- setNames(lapply(cohorts, function(cohort) {
  t(scaled[[cohort]][signature, , drop = FALSE])
}), cohorts)

for (cohort in cohorts) {
  if (!identical(rownames(m2efm_x[[cohort]]), rownames(clinical[[cohort]]))) {
    stop(cohort, ": model matrix/clinical order mismatch.")
  }
}

cat("Fitting final full-TCGA M2EFM molecular model...\n")
tcga_y <- survival::Surv(
  clinical$TCGA[[time_col]],
  clinical$TCGA[[event_col]]
)

foldid <- sample(rep(1:10, length.out = nrow(m2efm_x$TCGA)))

molecular_fit <- glmnet::cv.glmnet(
  x = m2efm_x$TCGA,
  y = tcga_y,
  family = "cox",
  alpha = 0,
  foldid = foldid,
  standardize = FALSE,
  type.measure = "deviance",
  grouped = TRUE,
  cox.ties = "breslow"
)

molecular_risk <- setNames(lapply(cohorts, function(cohort) {
  as.numeric(predict(
    molecular_fit,
    newx = m2efm_x[[cohort]],
    s = "lambda.min",
    type = "response"
  ))
}), cohorts)

tcga_integrated <- data.frame(
  survival_outcome = tcga_y,
  molecular_risk = molecular_risk$TCGA,
  age.Dx = clinical$TCGA[[age_col]],
  pathologic_stage = factor(
    as.character(clinical$TCGA[[stage_col]]),
    levels = stage_levels
  ),
  check.names = FALSE
)

cat("Fitting final full-TCGA integrated Cox model...\n")
integrated_fit <- survival::coxph(
  survival_outcome ~ molecular_risk + age.Dx + pathologic_stage,
  data = tcga_integrated,
  ties = "breslow",
  x = TRUE,
  model = TRUE
)

final_risk <- setNames(lapply(cohorts, function(cohort) {
  newdata <- data.frame(
    molecular_risk = molecular_risk[[cohort]],
    age.Dx = clinical[[cohort]][[age_col]],
    pathologic_stage = factor(
      as.character(clinical[[cohort]][[stage_col]]),
      levels = stage_levels
    ),
    check.names = FALSE
  )
  as.numeric(predict(
    integrated_fit,
    newdata = newdata,
    type = "lp",
    reference = "sample"
  ))
}), cohorts)

risk_data_list <- list()
summary_list <- list()
cindex_list <- list()
hr_list <- list()
logrank_list <- list()
km_fits <- list()

for (cohort in cohorts) {
  groups <- assign_groups(final_risk[[cohort]])

  dat <- data.frame(
    sample_id = rownames(clinical[[cohort]]),
    cohort = cohort,
    time_days = clinical[[cohort]][[time_col]],
    event = clinical[[cohort]][[event_col]],
    age_Dx = clinical[[cohort]][[age_col]],
    pathologic_stage = as.character(clinical[[cohort]][[stage_col]]),
    molecular_risk = molecular_risk[[cohort]],
    final_risk_score = final_risk[[cohort]],
    q25_cutoff = groups$q25,
    q75_cutoff = groups$q75,
    risk_group = groups$group,
    stringsAsFactors = FALSE
  )
  dat$risk_group <- factor(dat$risk_group, levels = risk_levels)
  risk_data_list[[cohort]] <- dat

  summary_rows <- lapply(risk_levels, function(group_name) {
    idx <- dat$risk_group == group_name
    data.frame(
      cohort = cohort,
      risk_group = group_name,
      samples = sum(idx),
      events = sum(dat$event[idx]),
      median_risk_score = median(dat$final_risk_score[idx]),
      stringsAsFactors = FALSE
    )
  })
  summary_list[[cohort]] <- do.call(rbind, summary_rows)

  cindex_list[[cohort]] <- cindex_row(final_risk[[cohort]], clinical[[cohort]], cohort)

  km_fits[[cohort]] <- survival::survfit(
    survival::Surv(time_days, event) ~ risk_group,
    data = dat
  )

  group_cox <- hr_rows(dat, cohort)
  hr_list[[cohort]] <- group_cox$table
  logrank_list[[cohort]] <- logrank_row(dat, cohort)
}

risk_data <- do.call(rbind, risk_data_list)
rownames(risk_data) <- NULL
group_summary <- do.call(rbind, summary_list)
rownames(group_summary) <- NULL
cindex_table <- do.call(rbind, cindex_list)
rownames(cindex_table) <- NULL
hr_table <- do.call(rbind, hr_list)
rownames(hr_table) <- NULL
logrank_table <- do.call(rbind, logrank_list)
rownames(logrank_table) <- NULL

write.csv(
  risk_data,
  file.path(outdir, "13_figure3_risk_scores.csv"),
  row.names = FALSE
)
write.csv(
  group_summary,
  file.path(outdir, "13_figure3_risk_group_summary.csv"),
  row.names = FALSE
)
write.csv(
  cindex_table,
  file.path(outdir, "13_figure3_final_model_cindices.csv"),
  row.names = FALSE
)
write.csv(
  hr_table,
  file.path(outdir, "13_figure3_hazard_ratios.csv"),
  row.names = FALSE
)
write.csv(
  logrank_table,
  file.path(outdir, "13_figure3_logrank_tests.csv"),
  row.names = FALSE
)

molecular_coef <- as.matrix(coef(molecular_fit, s = "lambda.min"))
write.csv(
  data.frame(
    gene = rownames(molecular_coef),
    coefficient = as.numeric(molecular_coef[, 1]),
    stringsAsFactors = FALSE
  ),
  file.path(outdir, "13_figure3_molecular_coefficients.csv"),
  row.names = FALSE
)

integrated_summary <- summary(integrated_fit)
integrated_coef <- data.frame(
  term = rownames(integrated_summary$coefficients),
  coefficient = integrated_summary$coefficients[, "coef"],
  hazard_ratio = integrated_summary$conf.int[, "exp(coef)"],
  lower_95 = integrated_summary$conf.int[, "lower .95"],
  upper_95 = integrated_summary$conf.int[, "upper .95"],
  p_value = integrated_summary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)
write.csv(
  integrated_coef,
  file.path(outdir, "13_figure3_integrated_model_coefficients.csv"),
  row.names = FALSE
)

plot_figure <- function(file, device = c("png", "pdf")) {
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

  for (cohort in cohorts) {
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
      xlim = c(0, ten_year_days),
      ylim = c(0.25, 1.02),
      xlab = "Days",
      ylab = "Survival",
      xaxt = "n",
      yaxt = "n"
    )

    axis(
      1,
      at = c(0, 1000, 2000, 3000),
      labels = c("0", "1000", "2000", "3000"),
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

png_file <- file.path(outdir, "13_figure3_KM_M2EFM_ExpClin.png")
pdf_file <- file.path(outdir, "13_figure3_KM_M2EFM_ExpClin.pdf")
plot_figure(png_file, "png")
plot_figure(pdf_file, "pdf")

final_object <- list(
  molecular_fit = molecular_fit,
  integrated_fit = integrated_fit,
  signature_genes = signature,
  scale_parameters = params,
  foldid = foldid,
  molecular_risk = molecular_risk,
  final_risk = final_risk,
  risk_data = risk_data,
  group_summary = group_summary,
  cindices = cindex_table,
  hazard_ratios = hr_table,
  logrank_tests = logrank_table,
  settings = list(
    model = "M2EFM Exp+Clin",
    training = "full TCGA",
    seed = 1L,
    alpha = 0,
    standardize = FALSE,
    cox_ties = "breslow",
    lambda = "lambda.min",
    risk_groups = "cohort-specific Q25/Q75",
    plot_limit_days = ten_year_days
  )
)

saveRDS(
  final_object,
  file.path(outdir, "13_figure3_final_model.rds"),
  compress = "xz"
)

audit <- data.frame(
  metric = c(
    "input_file",
    "final_model",
    "training_cohort",
    "TCGA_samples",
    "Terunuma_samples",
    "Kao_samples",
    "common_genes",
    "M2EFM_genes",
    "lambda_min",
    "risk_group_rule",
    "plot_limit_days"
  ),
  value = c(
    normalizePath(input_file, winslash = "/", mustWork = TRUE),
    "M2EFM Exp+Clin",
    "full TCGA",
    nrow(clinical$TCGA),
    nrow(clinical$Terunuma),
    nrow(clinical$Kao),
    length(common_genes),
    length(signature),
    molecular_fit$lambda.min,
    "within each cohort: Low<Q25; Medium=Q25-Q75; High>Q75",
    ten_year_days
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(outdir, "13_figure3_audit.csv"),
  row.names = FALSE
)

capture.output(
  {
    cat("STEP 13 — FIGURE 3\n\n")
    cat("FINAL MODEL C-INDICES\n")
    print(cindex_table)
    cat("\nRISK GROUP SUMMARY\n")
    print(group_summary)
    cat("\nHAZARD RATIOS\n")
    print(hr_table)
    cat("\nLOG-RANK TESTS\n")
    print(logrank_table)
    cat("\nINTEGRATED MODEL COEFFICIENTS\n")
    print(integrated_coef)
    cat("\nAUDIT\n")
    print(audit)
    cat("\nSESSION INFO\n")
    sessionInfo()
  },
  file = file.path(outdir, "13_figure3_sessionInfo.txt")
)

cat("\n============================================================\n")
cat("STEP 13 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nFINAL MODEL C-INDICES\n")
print(cindex_table, row.names = FALSE)

cat("\nRISK GROUP SUMMARY\n")
print(group_summary, row.names = FALSE)

cat("\nHAZARD RATIOS\n")
print(hr_table, row.names = FALSE)

cat("\nLOG-RANK TESTS\n")
print(logrank_table, row.names = FALSE)

cat("\nFigure saved:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat("\nFinished successfully.\n")
