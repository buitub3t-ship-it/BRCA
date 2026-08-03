scale_features_reference <- function(feature_by_sample) {
  assert_packages("matrixStats")
  means <- rowMeans(feature_by_sample)
  ranges <- matrixStats::rowMaxs(feature_by_sample) -
    matrixStats::rowMins(feature_by_sample)
  ranges[!is.finite(ranges) | ranges == 0] <- 1
  scaled <- sweep(feature_by_sample, 1L, means, FUN = "-")
  scaled <- sweep(scaled, 1L, ranges, FUN = "/")
  list(values = scaled, means = means, ranges = ranges)
}

apply_feature_reference <- function(feature_by_sample, means, ranges) {
  missing_features <- setdiff(names(means), rownames(feature_by_sample))
  if (length(missing_features) > 0L) {
    stop("Validation data miss features: ",
         paste(head(missing_features, 20L), collapse = ", "))
  }
  feature_by_sample <- feature_by_sample[names(means), , drop = FALSE]
  scaled <- sweep(feature_by_sample, 1L, means, FUN = "-")
  sweep(scaled, 1L, ranges, FUN = "/")
}

prepare_model_data <- function(feature_by_sample, clinical) {
  ids <- intersect(colnames(feature_by_sample), rownames(clinical))
  clinical <- clinical[ids, , drop = FALSE]
  feature_by_sample <- feature_by_sample[, ids, drop = FALSE]
  keep <- complete_os_rows(clinical)
  clinical <- clinical[keep, , drop = FALSE]
  feature_by_sample <- feature_by_sample[, rownames(clinical), drop = FALSE]
  list(features = feature_by_sample, clinical = clinical)
}

fit_two_stage_m2efm <- function(
    sample_by_feature,
    clinical,
    nfolds = GLMNET_CV_FOLDS,
    seed = 1L) {

  assert_packages(c("glmnet", "survival"))
  if (!identical(rownames(sample_by_feature), rownames(clinical))) {
    stop("Training molecular and clinical rows are not aligned.")
  }

  outcome <- survival::Surv(
    clinical$OS_days,
    clinical$overall.survival.indicator
  )
  nfolds <- min(as.integer(nfolds), nrow(sample_by_feature))
  if (nfolds < 3L) stop("At least three CV folds are required.")

  set.seed(seed)
  initial_model <- glmnet::cv.glmnet(
    x = data.matrix(sample_by_feature),
    y = outcome,
    family = "cox",
    alpha = 0,
    nfolds = nfolds,
    standardize = FALSE
  )

  molecular_risk <- as.numeric(predict(
    initial_model,
    newx = data.matrix(sample_by_feature),
    s = "lambda.min",
    type = "response"
  ))

  final_data <- data.frame(
    OS_days = clinical$OS_days,
    event = clinical$overall.survival.indicator,
    pred = molecular_risk,
    pathologic_stage = factor(
      clinical$pathologic_stage,
      levels = STAGE_LEVELS
    ),
    age.Dx = clinical$age.Dx,
    row.names = rownames(clinical)
  )

  final_model <- survival::coxph(
    survival::Surv(OS_days, event) ~ pred + pathologic_stage + age.Dx,
    data = final_data,
    x = TRUE,
    model = TRUE
  )

  clinical_model <- survival::coxph(
    survival::Surv(OS_days, event) ~ pathologic_stage + age.Dx,
    data = final_data,
    x = TRUE,
    model = TRUE
  )

  structure(
    list(
      initial_model = initial_model,
      final_model = final_model,
      clinical_model = clinical_model,
      feature_names = colnames(sample_by_feature)
    ),
    class = "m2efm_os_fit"
  )
}

predict_two_stage_m2efm <- function(fit, sample_by_feature, clinical) {
  if (!all(fit$feature_names %in% colnames(sample_by_feature))) {
    missing <- setdiff(fit$feature_names, colnames(sample_by_feature))
    stop("Prediction matrix lacks: ", paste(head(missing, 20L), collapse = ", "))
  }
  sample_by_feature <- sample_by_feature[, fit$feature_names, drop = FALSE]
  if (!identical(rownames(sample_by_feature), rownames(clinical))) {
    stop("Prediction molecular and clinical rows are not aligned.")
  }

  molecular_risk <- as.numeric(predict(
    fit$initial_model,
    newx = data.matrix(sample_by_feature),
    s = "lambda.min",
    type = "response"
  ))

  newdata <- data.frame(
    pred = molecular_risk,
    pathologic_stage = factor(
      clinical$pathologic_stage,
      levels = STAGE_LEVELS
    ),
    age.Dx = clinical$age.Dx,
    row.names = rownames(clinical)
  )

  combined_risk <- as.numeric(predict(
    fit$final_model,
    newdata = newdata,
    type = "lp"
  ))
  clinical_risk <- as.numeric(predict(
    fit$clinical_model,
    newdata = newdata,
    type = "lp"
  ))

  data.frame(
    sample_id = rownames(clinical),
    molecular_risk = molecular_risk,
    combined_risk = combined_risk,
    clinical_risk = clinical_risk,
    stringsAsFactors = FALSE,
    row.names = rownames(clinical)
  )
}

harrell_c_index <- function(clinical, risk) {
  value <- survival::concordance(
    survival::Surv(
      clinical$OS_days,
      clinical$overall.survival.indicator
    ) ~ risk,
    reverse = TRUE
  )$concordance
  as.numeric(value)
}

monte_carlo_m2efm <- function(
    training_feature_by_sample,
    training_clinical,
    validation_sets = list(),
    n_splits = N_MONTE_CARLO_SPLITS,
    training_proportion = TRAINING_PROPORTION,
    base_seed = BASE_SEED,
    keep_models = FALSE) {

  prepared_train <- prepare_model_data(
    training_feature_by_sample,
    training_clinical
  )
  train_features <- prepared_train$features
  train_clin <- prepared_train$clinical

  prepared_validation <- lapply(validation_sets, function(v) {
    prepare_model_data(v$features, v$clinical)
  })

  n <- ncol(train_features)
  n_train <- floor(training_proportion * n)
  if (n_train < 20L || n - n_train < 10L) {
    stop("Training/test split is too small.")
  }

  metric_rows <- vector("list", n_splits * (1L + length(prepared_validation)))
  model_list <- if (isTRUE(keep_models)) vector("list", n_splits) else NULL
  row_counter <- 1L

  for (i in seq_len(n_splits)) {
    seed <- base_seed + i - 1L
    set.seed(seed)
    train_index <- sample(seq_len(n), size = n_train, replace = FALSE)
    test_index <- setdiff(seq_len(n), train_index)

    train_ids <- colnames(train_features)[train_index]
    test_ids <- colnames(train_features)[test_index]

    x_train <- t(train_features[, train_ids, drop = FALSE])
    x_test <- t(train_features[, test_ids, drop = FALSE])
    clin_train <- train_clin[train_ids, , drop = FALSE]
    clin_test <- train_clin[test_ids, , drop = FALSE]

    fit <- fit_two_stage_m2efm(
      x_train,
      clin_train,
      seed = seed
    )
    if (isTRUE(keep_models)) model_list[[i]] <- fit

    pred_test <- predict_two_stage_m2efm(fit, x_test, clin_test)
    metric_rows[[row_counter]] <- data.frame(
      split = i,
      cohort = "TCGA_test",
      model = c("M2EFM_Molecular", "M2EFM_ExpClin", "Clinical"),
      c_index = c(
        harrell_c_index(clin_test, pred_test$molecular_risk),
        harrell_c_index(clin_test, pred_test$combined_risk),
        harrell_c_index(clin_test, pred_test$clinical_risk)
      )
    )
    row_counter <- row_counter + 1L

    if (length(prepared_validation) > 0L) {
      for (name in names(prepared_validation)) {
        v <- prepared_validation[[name]]
        x_val <- t(v$features)
        clin_val <- v$clinical
        pred_val <- predict_two_stage_m2efm(fit, x_val, clin_val)

        metric_rows[[row_counter]] <- data.frame(
          split = i,
          cohort = name,
          model = c("M2EFM_Molecular", "M2EFM_ExpClin", "Clinical"),
          c_index = c(
            harrell_c_index(clin_val, pred_val$molecular_risk),
            harrell_c_index(clin_val, pred_val$combined_risk),
            harrell_c_index(clin_val, pred_val$clinical_risk)
          )
        )
        row_counter <- row_counter + 1L
      }
    }

    message("Completed split ", i, "/", n_splits)
  }

  metrics <- do.call(rbind, metric_rows[seq_len(row_counter - 1L)])
  rownames(metrics) <- NULL
  list(metrics = metrics, models = model_list)
}


summarize_cindex_metrics <- function(metrics) {
  required <- c("cohort", "model", "c_index")
  absent <- setdiff(required, names(metrics))
  if (length(absent) > 0L) {
    stop("Metrics table lacks: ", paste(absent, collapse = ", "))
  }

  keys <- unique(metrics[c("cohort", "model")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    cohort_i <- keys$cohort[i]
    model_i <- keys$model[i]
    x <- metrics$c_index[
      metrics$cohort == cohort_i & metrics$model == model_i
    ]
    data.frame(
      cohort = cohort_i,
      model = model_i,
      n = length(x),
      median = stats::median(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      q025 = unname(stats::quantile(x, 0.025, na.rm = TRUE)),
      q975 = unname(stats::quantile(x, 0.975, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$cohort, out$model), , drop = FALSE]
}

fit_final_and_validate <- function(
    training_feature_by_sample,
    training_clinical,
    validation_sets = list(),
    seed = 10001L) {

  train <- prepare_model_data(training_feature_by_sample, training_clinical)
  fit <- fit_two_stage_m2efm(
    t(train$features),
    train$clinical,
    seed = seed
  )

  predictions <- list()
  metrics <- list()

  pred_train <- predict_two_stage_m2efm(
    fit,
    t(train$features),
    train$clinical
  )
  pred_train$OS_days <- train$clinical$OS_days
  pred_train$event <- train$clinical$overall.survival.indicator
  predictions$TCGA_train <- pred_train
  metrics$TCGA_train <- data.frame(
    cohort = "TCGA_train",
    model = c("M2EFM_Molecular", "M2EFM_ExpClin", "Clinical"),
    c_index = c(
      harrell_c_index(train$clinical, pred_train$molecular_risk),
      harrell_c_index(train$clinical, pred_train$combined_risk),
      harrell_c_index(train$clinical, pred_train$clinical_risk)
    )
  )

  for (name in names(validation_sets)) {
    v <- prepare_model_data(
      validation_sets[[name]]$features,
      validation_sets[[name]]$clinical
    )
    pred <- predict_two_stage_m2efm(fit, t(v$features), v$clinical)
    pred$OS_days <- v$clinical$OS_days
    pred$event <- v$clinical$overall.survival.indicator
    predictions[[name]] <- pred
    metrics[[name]] <- data.frame(
      cohort = name,
      model = c("M2EFM_Molecular", "M2EFM_ExpClin", "Clinical"),
      c_index = c(
        harrell_c_index(v$clinical, pred$molecular_risk),
        harrell_c_index(v$clinical, pred$combined_risk),
        harrell_c_index(v$clinical, pred$clinical_risk)
      )
    )
  }

  list(
    fit = fit,
    metrics = do.call(rbind, metrics),
    predictions = predictions
  )
}

save_metric_boxplot <- function(metrics, path, title) {
  grDevices::pdf(path, width = 10, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::boxplot(
    c_index ~ interaction(cohort, model, sep = "\n"),
    data = metrics,
    las = 2,
    ylab = "Harrell C-index",
    xlab = "",
    main = title,
    ylim = c(0.45, 0.95)
  )
  graphics::abline(h = 0.70, lty = 2)
}
