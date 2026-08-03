# Plotting helpers for author-style OS-M2EFM figures.
# These functions use only base graphics and the survival package.

open_plot_device <- function(path, width, height, res = 300L) {
  ext <- tolower(tools::file_ext(path))
  if (identical(ext, "pdf")) {
    grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  } else if (identical(ext, "png")) {
    grDevices::png(
      path,
      width = round(width * res),
      height = round(height * res),
      res = res,
      type = if (.Platform$OS.type == "windows") "windows" else "cairo"
    )
  } else {
    stop("Unsupported plot extension: ", ext, ". Use PDF or PNG.")
  }
}

plot_model_label <- function(model, context = c("internal", "external")) {
  context <- match.arg(context)
  if (context == "internal") {
    labels <- c(
      Clinical = "Clinical",
      M2EFM_Molecular = "M2EFM Meth+Exp",
      M2EFM_MethExp = "M2EFM Meth+Exp",
      M2EFM_ExpClin = "M2EFM Meth+Exp+Clin",
      M2EFM_MethExpClin = "M2EFM Meth+Exp+Clin"
    )
  } else {
    labels <- c(
      Clinical = "Clinical",
      M2EFM_Molecular = "M2EFM Exp",
      M2EFM_Exp = "M2EFM Exp",
      M2EFM_ExpClin = "M2EFM Exp+Clin"
    )
  }
  out <- unname(labels[model])
  out[is.na(out)] <- model[is.na(out)]
  out
}

save_author_figure2a <- function(metrics, path) {
  dat <- metrics[metrics$cohort == "TCGA_test", , drop = FALSE]
  dat$display_model <- plot_model_label(dat$model, "internal")
  model_order <- c("Clinical", "M2EFM Meth+Exp", "M2EFM Meth+Exp+Clin")
  dat <- dat[dat$display_model %in% model_order, , drop = FALSE]
  dat$display_model <- factor(dat$display_model, levels = model_order)

  groups <- lapply(model_order, function(x) {
    dat$c_index[dat$display_model == x]
  })

  open_plot_device(path, width = 7.2, height = 5.6)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(7.2, 4.5, 2.8, 1.2), las = 1)

  graphics::boxplot(
    groups,
    names = model_order,
    col = c("grey82", "grey65", "grey35"),
    border = "black",
    outline = FALSE,
    ylim = c(0.45, 0.95),
    ylab = "C-index",
    xlab = "",
    main = "A) TCGA internal validation",
    las = 2,
    cex.axis = 0.90
  )
  graphics::abline(h = 0.70, lty = 2, col = "grey40")
  medians <- vapply(groups, stats::median, numeric(1), na.rm = TRUE)
  graphics::points(seq_along(medians), medians, pch = 19, cex = 0.75)
  graphics::mtext("100 random 70/30 splits", side = 3, line = 0.2, cex = 0.82)
}

save_author_figure2b <- function(metrics, path) {
  cohorts <- c("TCGA_test", "Terunuma", "Kao")
  cohort_labels <- c(TCGA_test = "TCGA", Terunuma = "Terunuma", Kao = "Kao")
  models <- c("Clinical", "M2EFM_Molecular", "M2EFM_ExpClin")
  model_labels <- c("Clinical", "M2EFM Exp", "M2EFM Exp+Clin")
  model_cols <- c("grey82", "grey60", "grey30")

  box_data <- list()
  positions <- numeric()
  box_cols <- character()
  counter <- 1L
  for (i in seq_along(cohorts)) {
    for (j in seq_along(models)) {
      box_data[[counter]] <- metrics$c_index[
        metrics$cohort == cohorts[i] & metrics$model == models[j]
      ]
      positions[counter] <- (i - 1L) * 4L + j
      box_cols[counter] <- model_cols[j]
      counter <- counter + 1L
    }
  }

  open_plot_device(path, width = 9.2, height = 5.6)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4.8, 4.5, 3.0, 1.2), las = 1)

  graphics::boxplot(
    box_data,
    at = positions,
    names = rep("", length(box_data)),
    col = box_cols,
    border = "black",
    outline = FALSE,
    ylim = c(0.45, 0.95),
    xlim = c(0.3, 11.7),
    ylab = "C-index",
    xlab = "",
    main = "B) Expression-only internal and external validation",
    xaxt = "n"
  )
  graphics::axis(
    side = 1,
    at = c(2, 6, 10),
    labels = unname(cohort_labels[cohorts]),
    tick = FALSE,
    line = 0.4
  )
  graphics::abline(h = 0.70, lty = 2, col = "grey40")
  graphics::abline(v = c(4, 8), lty = 3, col = "grey75")
  graphics::legend(
    "topright",
    legend = model_labels,
    fill = model_cols,
    border = "black",
    bty = "n",
    cex = 0.82
  )
  graphics::mtext("100 TCGA random 70/30 splits; fixed external cohorts", side = 3, line = 0.2, cex = 0.82)
}

assign_author_risk_groups <- function(prediction_table) {
  required <- c("sample_id", "combined_risk", "OS_days", "event")
  missing <- setdiff(required, names(prediction_table))
  if (length(missing) > 0L) {
    stop("Prediction table lacks: ", paste(missing, collapse = ", "))
  }

  q <- stats::quantile(
    prediction_table$combined_risk,
    probs = c(0.25, 0.75),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
  group <- ifelse(
    prediction_table$combined_risk < q[1],
    "Low",
    ifelse(prediction_table$combined_risk > q[2], "High", "Medium")
  )

  out <- prediction_table
  out$risk_group <- factor(group, levels = c("Low", "Medium", "High"))
  out$risk_q25 <- q[1]
  out$risk_q75 <- q[2]
  out
}

risk_group_logrank_p <- function(data) {
  fit <- survival::survdiff(
    survival::Surv(OS_days, event) ~ risk_group,
    data = data
  )
  stats::pchisq(fit$chisq, df = length(fit$n) - 1L, lower.tail = FALSE)
}

risk_group_hr_table <- function(data, cohort) {
  fit <- survival::coxph(
    survival::Surv(OS_days, event) ~ risk_group,
    data = data
  )
  sm <- summary(fit)
  if (nrow(sm$coefficients) == 0L) {
    return(data.frame())
  }
  data.frame(
    cohort = cohort,
    comparison = sub("risk_group", "", rownames(sm$coefficients), fixed = TRUE),
    reference = "Low",
    HR = unname(sm$coefficients[, "exp(coef)"]),
    CI_low = unname(sm$conf.int[, "lower .95"]),
    CI_high = unname(sm$conf.int[, "upper .95"]),
    p_value = unname(sm$coefficients[, "Pr(>|z|)"]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

save_author_figure3 <- function(predictions, path) {
  required_cohorts <- c("TCGA_train", "Terunuma", "Kao")
  absent <- setdiff(required_cohorts, names(predictions))
  if (length(absent) > 0L) {
    stop("Final predictions lack: ", paste(absent, collapse = ", "))
  }

  cohort_titles <- c(TCGA_train = "A) TCGA", Terunuma = "B) Terunuma", Kao = "C) Kao")
  risk_cols <- c(Low = "#2B8C3E", Medium = "#2C7FB8", High = "#D7301F")
  lty <- c(Low = 1, Medium = 2, High = 3)

  grouped <- lapply(predictions[required_cohorts], assign_author_risk_groups)

  open_plot_device(path, width = 12.0, height = 4.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(4.3, 4.2, 3.1, 1.0), las = 1)

  for (cohort in required_cohorts) {
    dat <- grouped[[cohort]]
    dat$OS_years <- dat$OS_days / 365.25
    fit <- survival::survfit(
      survival::Surv(OS_years, event) ~ risk_group,
      data = dat
    )

    graphics::plot(
      fit,
      col = unname(risk_cols),
      lty = unname(lty),
      lwd = 2.2,
      mark.time = TRUE,
      conf.int = FALSE,
      xlim = c(0, 10),
      ylim = c(0, 1),
      xlab = "Years",
      ylab = "Survival probability",
      main = cohort_titles[cohort],
      yaxs = "i"
    )
    p <- risk_group_logrank_p(dat)
    graphics::legend(
      "bottomleft",
      legend = levels(dat$risk_group),
      col = unname(risk_cols),
      lty = unname(lty),
      lwd = 2.2,
      bty = "n",
      cex = 0.82
    )
    graphics::text(
      x = 9.7,
      y = 0.96,
      labels = paste0("Log-rank p = ", format.pval(p, digits = 2, eps = 1e-4)),
      adj = c(1, 1),
      cex = 0.78
    )
  }

  invisible(grouped)
}

baseline_hazard_at <- function(cox_model, horizon_days) {
  bh <- survival::basehaz(cox_model, centered = TRUE)
  eligible <- which(bh$time <= horizon_days)
  if (length(eligible) == 0L) return(0)
  bh$hazard[max(eligible)]
}

rank_calibration_bins <- function(x, n_bins = 5L) {
  valid <- which(is.finite(x))
  out <- rep(NA_integer_, length(x))
  if (length(valid) == 0L) return(out)
  ordered <- valid[order(x[valid], seq_along(valid))]
  out[ordered] <- pmin(
    n_bins,
    ceiling(seq_along(ordered) * n_bins / length(ordered))
  )
  out
}

observed_event_probability <- function(OS_days, event, horizon_days) {
  sf <- survival::survfit(survival::Surv(OS_days, event) ~ 1)
  surv_h <- summary(sf, times = horizon_days, extend = TRUE)$surv
  if (length(surv_h) == 0L || !is.finite(surv_h[1])) return(NA_real_)
  1 - surv_h[1]
}

build_calibration_table <- function(
    final_model,
    predictions,
    horizon_days = 5 * 365.25,
    n_bins = 5L) {

  h0 <- baseline_hazard_at(final_model, horizon_days)
  cohort_map <- c(TCGA_train = "TCGA", Kao = "Kao", Terunuma = "Terunuma")
  rows <- list()
  counter <- 1L

  for (key in names(cohort_map)) {
    dat <- predictions[[key]]
    if (is.null(dat)) next
    dat$predicted_event <- 1 - exp(-h0 * exp(dat$combined_risk))
    dat$calibration_bin <- rank_calibration_bins(dat$predicted_event, n_bins)

    for (bin in sort(unique(dat$calibration_bin[!is.na(dat$calibration_bin)]))) {
      d <- dat[dat$calibration_bin == bin, , drop = FALSE]
      rows[[counter]] <- data.frame(
        cohort = unname(cohort_map[key]),
        bin = bin,
        n = nrow(d),
        mean_predicted_event = mean(d$predicted_event, na.rm = TRUE),
        observed_event = observed_event_probability(
          d$OS_days,
          d$event,
          horizon_days
        ),
        horizon_years = horizon_days / 365.25,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }

  do.call(rbind, rows)
}

save_author_figureS1 <- function(calibration_table, path) {
  cohorts <- c("TCGA", "Kao", "Terunuma")
  cohort_cols <- c(TCGA = "black", Kao = "#2C7FB8", Terunuma = "#D7301F")
  cohort_pch <- c(TCGA = 16, Kao = 17, Terunuma = 15)
  max_value <- max(
    calibration_table$mean_predicted_event,
    calibration_table$observed_event,
    na.rm = TRUE
  )
  axis_max <- min(1, max(0.75, ceiling(max_value * 20) / 20))

  open_plot_device(path, width = 6.2, height = 6.0)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4.8, 4.8, 3.2, 1.2), las = 1)

  graphics::plot(
    NA_real_,
    NA_real_,
    xlim = c(0, axis_max),
    ylim = c(0, axis_max),
    xlab = "Predicted event probability",
    ylab = "Observed event frequency",
    main = "M2EFM Exp+Clin calibration at 5 years",
    xaxs = "i",
    yaxs = "i"
  )
  graphics::abline(0, 1, lty = 2, col = "grey50")

  for (cohort in cohorts) {
    d <- calibration_table[calibration_table$cohort == cohort, , drop = FALSE]
    d <- d[order(d$mean_predicted_event), , drop = FALSE]
    graphics::lines(
      d$mean_predicted_event,
      d$observed_event,
      type = "b",
      lwd = 2,
      pch = cohort_pch[cohort],
      col = cohort_cols[cohort]
    )
  }
  graphics::legend(
    "topleft",
    legend = cohorts,
    col = unname(cohort_cols[cohorts]),
    pch = unname(cohort_pch[cohorts]),
    lwd = 2,
    bty = "n"
  )
  graphics::mtext("Grouped Kaplan-Meier calibration", side = 3, line = 0.2, cex = 0.82)
}
