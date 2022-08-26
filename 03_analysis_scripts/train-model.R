#' Script to train, evaluate, and make predictions of extinction risk

#' EXPECTED INPUTS:
#'  - `method`: the name of a modelling method to use, e.g `bart`
#'  - `output_dir`: path to a directory to save outputs to
#'  - `model_dir`: path to a directory to save trained models to
#'  - `method_dir`: path to a directory containing scripts that specify each method
#'  - `predictor_file`: path to a file containing predictors calculated for a set of species
#' 
#' EXAMPLE CLI:
#'  Rscript 03_analysis_scripts/train-model.R --predictor_file=01_raw_data/angio_pred_v1.csv --method=bart
#' 
#' EXAMPLE SOURCE:
#'  output_dir <- "output"
#'  model_dir <- "output"
#'  method_dir <- "04_model_definitions"
#'  predictor_file <- "01_raw_data/angio_pred_v1.csv"
#'  random_seed <- 1989
#'  method <- "bart"
#'  
#'  source("03_analysis_scripts/train-model.R")
#' 

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(tidyverse)   # packages for data handling
shhlibrary(tidymodels)  # packages for model training and evaluation
shhlibrary(ROCR)        # evaluate classification thresholds
shhlibrary(multidplyr)  # run operations on chunks of a data frame in parallel
shhlibrary(doParallel)  # set up parallel processing
shhlibrary(git2r)       # run git commands in R
shhlibrary(cli)         # nice formatting for CLI

source("R/model-functions.R")

# CLI ----
cli_h1("Running extinction risk prediction method")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="output",
    method_dir="04_model_definitions",
    model_dir="output",
    random_seed=1989
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  method <- args$method
  predictor_file <- args$predictor_file
  
  random_seed <- args$random_seed
  
  output_dir <- args$output_dir
  method_dir <- args$method_dir
  model_dir <- args$model_dir
}

if (! exists("random_seed")) {
  random_seed <- NULL
}

if (! exists("predictor_file", mode="character")) {
  cli_abort(c(
    "no path to predictors provided",
    "x"="You must provide a path to a file of species-level predictors as {.var predictor_file}."
  ))
}

if (! exists("method", mode="character")) {
  cli_abort(c(
    "no method provided",
    "x"="You must specify which automated assessment method to evaluate as {.var method}."
  ))
}

if (! exists("output_dir")) {
  cli_abort(c(
    "no path to save results provided",
    "x"="You must provide the save path as the variable {.var output_dir}."
  ))
}

if (! exists("model_dir")) {
  cli_abort(c(
    "no path to save trained models provided",
    "x"="You must provide the save path as the variable {.var model_dir}."
  ))
}

available_methods <- 
  method_dir |>
  list.files() |>
  str_remove("\\.R")

if (! method %in% available_methods) {
  cli_abort(c(
    "unrecognised method",
    x="No implementation for method found, implement it yourself or use one of {.var {available_methods}}."
  ))
}

source(file.path(method_dir, paste0(method, ".R")))

dir.create(output_dir, showWarnings=FALSE)
dir.create(model_dir, showWarnings=FALSE)

latest_hash <- revparse_single(".", "HEAD")$sha
now <- format(Sys.time(), "%Y%m%d-%H%M%S")
name <- paste(method, latest_hash, now, sep="-")

cli_alert_info("Evaluating {.strong {method}} method on {.file {predictor_file}}")
cli_alert_info("Saving results to {.file {output_dir}}")

# load predictors ----
predictors <- read_csv(predictor_file, show_col_types=FALSE, progress=FALSE)

# prepare predictors ----
predictors <- filter(predictors, ! category %in% c("EW", "EX"))

old_codes_map <- c("LR/cd"="NT", "LR/lc"="LC", "LR/nt"="NT")
predictors$category <- recode(predictors$category, !!! old_codes_map)

predictors$humphreys_lifeform <- factor(predictors$humphreys_lifeform)
predictors$climate_description <- factor(predictors$climate_description)

labelled <- filter(predictors, ! is.na(category), category != "DD")
unlabelled <- filter(predictors, is.na(category) | category == "DD")

labelled$obs <- ifelse(labelled$category %in% c("LC", "NT"), "not threatened", "threatened")
labelled$obs <- factor(labelled$obs, levels=c("not threatened", "threatened"))

unlabelled$obs <- factor(NA, levels=levels(labelled$obs), ordered=TRUE)

pct_threat <- mean(labelled$obs == 'threatened')

cli_alert_info("Evaluating the method on {.strong {nrow(labelled)}} examples ({.strong {label_percent()(pct_threat)}} threatened).")

# define data budget ----
set.seed(random_seed)

random_cv <- nested_cv(labelled, inside=vfold_cv(v=3), outside=vfold_cv(v=5))
family_cv <- nested_cv(labelled, inside=vfold_cv(v=3), 
                       outside=group_vfold_cv(group=family, v=5))

# define evaluation metrics ----
# hyperparameter tuning (inner folds)
tune_metrics <- metric_set(roc_auc, mn_log_loss)

# model evaluation (outer folds)
eval_metrics <- metric_set(accuracy, sensitivity, specificity, j_index)

# set up cluster ----
ncores <- parallelly::availableCores()
cli_alert_info("Using {.strong {ncores}} cores")

# this sets up a cluster for hyperparameter tuning
cl <- makePSOCKcluster(ncores)
registerDoParallel(cl)

# this sets up a cluster for mapping chunks of data frames (model evaluation)
cluster <- new_cluster(ncores)
cluster_library(cluster, packages=c("dplyr", "tidymodels", "ROCR"))

# set up model ----
model_spec <- specify_model()
model_recipe <- specify_recipe(labelled)

wf <- 
  workflow() |>
  add_model(model_spec) |>
  add_recipe(model_recipe)

# tune hyperparameters if needed ----
if (exists("hparam_grid")) {
  cli_alert_info("Tuning on inner resamples to find best hyperparameters")
  
  random_cv <- 
    random_cv |>
    tune_hyperparameters(wf, grid=hparam_grid, metrics=tune_metrics, parallel=T)
  
  family_cv <- 
    family_cv |>
    tune_hyperparameters(wf, grid=hparam_grid, metrics=tune_metrics, parallel=T)
  
} else {
  cli_alert_info("No hyperparameters to tune, just evaluating on outer folds")
  random_cv <- 
    random_cv |>
    mutate(.workflow=list(wf))
  
  family_cv <- 
    family_cv |>
    mutate(.workflow=list(wf))
}

# evaluate model by cv ----
random_results <- evaluate_model(random_cv, metrics=eval_metrics, cluster=cluster)
cli_alert_success("Evaluated model using random CV on assessed species")

family_results <- evaluate_model(family_cv, metrics=eval_metrics, cluster=cluster)
cli_alert_success("Evaluated model using taxonomic-block CV on assessed species")

# save results
random_performance <- extract_performance(random_results)
family_performance <- extract_performance(family_results)

write_csv(random_performance, file.path(output_dir, paste0(name, "-random-performance.csv")))
write_csv(family_performance, file.path(output_dir, paste0(name, "-family-performance.csv")))

random_test_pred <- extract_predictions(random_results)
family_test_pred <- extract_predictions(family_results)

write_csv(random_test_pred, file.path(output_dir, paste0(name, "-random-test-preds.csv")))
write_csv(random_test_pred, file.path(output_dir, paste0(name, "-family-test-preds.csv")))

# tune and fit final model ----
final_wf <- wf

if (exists("hparam_grid")) {
  cli_alert_info("Tuning hyperparameters of final model")
  
  final_tune <- tune_grid(
    final_wf, 
    vfold_cv(labelled, v=5), 
    grid=hparam_grid, 
    metrics=tune_metrics,
    control=control_grid(parallel_over="everything")
  )
  best_params <- select_best(final_tune, metric="roc_auc")
  final_wf <- finalize_workflow(final_wf, best_params)
}

cli_alert_info("Training final model on all labelled data")
fit_wf <- fit(final_wf, labelled)

if (method == "bart") {
  # have to 'touch' the trees to be able to save them for predictions
  invisible(fit_wf$fit$fit$fit$fit$state)  
}


write_rds(fit_wf, file.path(model_dir, paste0(name, ".rds")))

# make predictions ----
cli_alert_info("Generating predictions for all species")

if (model_spec$engine == "dbarts") {
  predictions <- predict_bart(fit_wf, labelled, unlabelled)
  
  write_csv(predictions$ppd, file.path(output_dir, paste0(name, "-ppd-samples.csv")))
  write_csv(predictions$ev, file.path(output_dir, paste0(name, "-ev-samples.csv")))
}

predictions <- predict_classes(fit_wf, labelled, unlabelled)
write_csv(predictions, file.path(output_dir, paste0(name, "-final-predictions.csv")))

cli_alert_success("Finished training, evaluations, predictions, and saving outputs!")
