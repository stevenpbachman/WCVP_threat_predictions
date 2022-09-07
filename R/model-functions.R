#' Functions to help train and evaluate the model

#' Tune the hyperparameters of a model across nested CV folds.
#' 
#' @param splits an `rsample` object holding nested CV folds
#' @param wf an unfinalised `workflow` object with parameters to tune
#' @param grid a data frame specifying parameter values to try
#' @param metrics a metric or set of metrics from to evaluate parameter values
#' @param parallel whether or not to use parallel processing, you must have set
#'  up a cluster already to use this
#'  
#' @return the original `rsample` object modified with finalised workflows using the
#'  best parameter values from each split
#'  
tune_hyperparameters <- function(splits, wf, grid, metrics, parallel=FALSE) {
  splits <- 
    splits |>
    rowwise() |>
    mutate(tuning=list(tune_over_folds(inner_resamples, wf, grid, metrics,
                                       parallel=parallel)))
  
  splits |>
    rowwise() |>
    mutate(best=list(select_best(tuning, metric="roc_auc"))) |>
    mutate(.workflow=list(finalize_workflow(wf, best))) |>
    select(-tuning) |>
    unnest(cols=c(best))
}

#' Evaluate the performance of a workflow for different hyperparameter values.
#' 
#' @param folds an `rsplit` object.
#' @param wf an unfinalised `workflow` object with parameters to tune
#' @param grid a data frame specifying parameter values to try
#' @param metrics a metric or set of metrics from to evaluate parameter values
#' @param parallel whether or not to use parallel processing, you must have set
#'  up a cluster already to use this
#' 
#' @return An updated version of the `folds` object with extra columns for the 
#'  tuning results.
#'
tune_over_folds <- function(folds, wf, grid, metrics, parallel=FALSE) {
  if (parallel) {
    control <- control_grid(parallel_over="everything")
  } else {
    control <- control_grid()
  }
  
  tune_grid(
    wf, 
    folds, 
    grid=hparam_grid, 
    metrics=tune_metrics,
    control=control
  )
}

#' Fit and evaluate the performance of a model over cross-validation folds.
#'
#' @param splits an `rsample::rsets` object with a `.workflow` column
#'  holding the finalised workflows. These will be the same if no hyperparameter
#'  tuning on nested resamples has been performed.
#' @param metrics a metric or set of metrics to measure the performance of the model.
#' @param cluster a cluster to carry out parallel processing using `multidplyr`
evaluate_model <- function(splits, metrics, cluster=NULL) {
  if (! is.null(cluster)) {
    cluster_assign(
      cluster,
      #functions
      metrics=metrics,
      last_fit_threshold=last_fit_threshold,
      choose_threshold=choose_threshold,
      permutation_importance=permutation_importance
    )  
  }
  
  results <-
    splits |>
    rowwise()
  
  if (! is.null(cluster)) {
    results <- partition(results, cluster)
  }
  
  results <- 
    results |>
    mutate(.fit=list(last_fit_threshold(.workflow, splits, metrics=metrics))) |>
    mutate(.importance=list(permutation_importance(.fit$.fit, splits, .threshold=.threshold)))
  
  if (! is.null(cluster)) {
    results <- collect(results)
  }
  
  results
}

#' Calculate predictor importance by random permutations.
#'
#' @param object A fitted `workflow` object.
#' @param split A single CV fold from `rsample`.
#' 
#' @return A data frame of the predictor importance measured as the mean decrease
#'  in accuracy after randomly permuting each predictor in turn.
#'  
permutation_importance <- function(object, split, .threshold=0.5) {
  
  trained_rec <- extract_recipe(object)
  fit_obj <- extract_fit_engine(object)
  
  newdata <- bake(trained_rec, assessment(split))
  if (class(fit_obj) == "bart") {
    fcn <- function(object, newdata) {
      p <- colMeans(predict(object, newdata))
      ifelse(p > .threshold, "threatened", "not threatened")
    }
  } else {
    fcn <- function(object, newdata) {
      p <- predict(object, newdata, type="prob")
      ifelse(p > .threshold, "threatened", "not threatened")
    }
  }
  
  vip::vi_permute(
    fit_obj,
    train=newdata,
    target="obs",
    metric="accuracy",
    pred_wrapper=fcn,
    nsim=50
  )
}

#' Fit and evaluate a model with additional threshold tuning.
#' 
#' Fits the model on the training set of a single CV fold, picks the classification
#' threshold that maximises the TSS, and then evaluates the model performance with 
#' that threshold using class-based performance metrics.
#' 
#' @param wf A finalised `workflow` object.
#' @param split A single CV fold or train/test split as an `rsample::rset` object.
#' @param metrics A single or set of class-based classification metrics, e.g. accuracy.
#' 
last_fit_threshold <- function(wf, split, metrics) {
  .fit <- fit(wf, analysis(split))
  .pred <- augment(.fit, assessment(split))
  
  .threshold <- choose_threshold(.pred$.pred_threatened, .pred$obs)
  .pred$.pred_class <- ifelse(.pred$.pred_threatened > .threshold$best,
                              "threatened", "not threatened")
  .pred$.pred_class <- factor(.pred$.pred_class, levels=levels(.pred$obs))
  
  .metrics <- metrics(.pred, truth=obs, estimate=.pred_class, event_level="second")
  
  list(
    .fit=.fit,
    .predictions=.pred,
    .metrics=.metrics,
    .threshold=.threshold$best
  )
}

#' Given predicted class probabilites, choose the classification threshold that maximises
#' the TSS.
#' 
#' @preds A vector of predicted class probabilities.
#' @labels A vector of observed class labels corresponding to the predictions.
#' 
#' @return A list holding the threshold values tried (`alpha`), the TSS, and the 
#'  chosen threshold (`best`).
#'  
choose_threshold <- function(preds, labels) {
  rocr_preds <- prediction(preds, labels, label.order=c("not threatened", "threatened"))
  rocr_perf <- performance(rocr_preds, "sens", "spec")
  
  tss <- rocr_perf@x.values[[1]] + rocr_perf@y.values[[1]] - 1
  alpha <- rocr_perf@alpha.values[[1]]
  list(
    alpha=alpha,
    tss=tss,
    best=alpha[which.max(tss)]
  )
}

#' Extract the cross-validated performance of a model from a nested data frame.
#' 
#' @param results A data frame of modelling results after running `evaluate_model`.
#' 
extract_performance <- function(results) {
  results |>
    rowwise() |>
    mutate(.metrics=list(.fit$.metrics)) |>
    mutate(.threshold=.fit$.threshold) |>
    select(id, .metrics, .threshold) |>
    unnest(.metrics)
}

#' Extract the cross-validated test-set predictions from a nested data frame.
#' 
#' @param results A data frame of modelling results after running `evaluate_model`.
#' 
extract_predictions <- function(results) {
  results |>
    rowwise() |>
    mutate(.preds=list(.fit$.predictions)) |>
    select(id, .preds) |>
    unnest(.preds)
}

#' Extract the cross-validated permutation importance from a nested data frame.
#' 
#' @param results A data frame of modelling results after running `evaluate_model`.
#' 
extract_importance <- function(results) {
  results |>
    select(id, .importance) |>
    unnest(.importance) |>
    select(variable=Variable, mean_decrease_accuracy=Importance)
}


#' Sample posterior predictions for labelled and unlabelled data using a BART model.
#'
#' @param wf a fit `workflow` object
#' @param labelled a data frame of labelled data for predictions
#' @param unlabelled a data frame of unlabelled data for predictions
predict_bart <- function(wf, labelled, unlabelled) {
  # extracting fit object and recipe so we can sample posterior
  trained_rec <- extract_recipe(wf)
  fit_obj <- extract_fit_engine(wf)
  
  # process data for predictions
  labelled_processed <- bake(trained_rec, new_data=labelled)
  unlabelled_processed <- bake(trained_rec, new_data=unlabelled)
  
  # sample posterior probability of threat
  pred_labelled <- predict(fit_obj, newdata=labelled_processed, type="ev")
  pred_unlabelled <- predict(fit_obj, newdata=unlabelled_processed, type="ev")
  
  ev_samples <-
    pred_labelled |>
    t() |>
    as_tibble() |>
    mutate(set="labelled", plant_name_id=labelled$plant_name_id) |>
    bind_rows(
      pred_unlabelled |>
        t() |>
        as_tibble() |>
        mutate(set="unlabelled", plant_name_id=unlabelled$plant_name_id)
    )
  
  # sample posterior predictive (so we can get estimated proportion threatened)
  ppd_labelled <- predict(fit_obj, newdata=labelled_processed, type="ppd")
  ppd_unlabelled <- predict(fit_obj, newdata=unlabelled_processed, type="ppd")
  
  ppd_samples <-
    ppd_labelled |>
    t() |>
    as_tibble() |>
    mutate(set="labelled", plant_name_id=labelled$plant_name_id) |>
    bind_rows(
      ppd_unlabelled |> 
        t() |>
        as_tibble() |>
        mutate(set="unlabelled", plant_name_id=unlabelled$plant_name_id)
    )
  
  list(
    ev=ev_samples,
    ppd=ppd_samples
  )
}

#' Make predictions for labelled and unlabelled data.
#'
#' @param wf a fit `workflow` object
#' @param labelled a data frame of labelled data for predictions
#' @param unlabelled a data frame of unlabelled data for predictions
predict_classes <- function(wf, labelled, unlabelled, threshold=NULL) {
  labelled_pred <- augment(wf, new_data=labelled)
  if (is.null(threshold)) {
    threshold <- choose_threshold(labelled_pred$.pred_threatened, labelled_pred$obs)
    threshold <- threshold$best
  }
  
  unlabelled_pred <- augment(wf, new_data=unlabelled)
  
  labelled_pred |>
    mutate(set="labelled") |>
    bind_rows(
      unlabelled_pred |>
        mutate(set="unlabelled")
    ) |>
    mutate(.pred_class=ifelse(.pred_threatened > threshold,
                              "threatened", "not threatened"),
           .threshold=threshold) |>
    mutate(.pred_class=factor(.pred_class, levels=levels(obs)))
}
