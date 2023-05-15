#' Set up the hyperparameters, cross-validation folds, and pre-processor for 
#' training a model.


#' EXPECTED INPUTS:
#'  - `method`: the name of a modelling method to use, e.g `bart`
#'  - `output_dir`: path to a directory to save outputs to
#'  - `model_dir`: path to a directory to save trained models to
#'  - `method_dir`: path to a directory containing scripts that specify each method
#'  - `predictor_file`: path to a file containing predictors calculated for a set of species 
#'
#' EXAMPLE CLI:
#'  Rscript 03_analysis_scripts/setup-model.R --predictor_file=01_raw_data/angio_pred_v1.csv --method=bart
#' 
#' EXAMPLE SOURCE:
#'  output_dir <- "output"
#'  method_dir <- "04_model_definitions"
#'  predictor_file <- "01_raw_data/angio_pred_v1.csv"
#'  random_seed <- 1989
#'  method <- "bart"
#'  
#'  source("03_analysis_scripts/train-model.R")
#' 

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(dbarts)      # for a bart model if needed
shhlibrary(tidyverse)   # packages for data handling
shhlibrary(tidymodels)  # packages for model training and evaluation
shhlibrary(ROCR)        # evaluate classification thresholds
shhlibrary(multidplyr)  # run operations on chunks of a data frame in parallel
shhlibrary(doParallel)  # set up parallel processing
shhlibrary(git2r)       # run git commands in R
shhlibrary(cli)         # nice formatting for CLI

source("R/utility-functions.R")
source("R/model-functions.R")
# CLI ----
cli_h1("Setting up an extinction prediction model training run")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="05_outputs",
    method_dir="04_model_definitions",
    mode="eval",
    threshold=NULL,
    random_seed=1989,
    force_commits=TRUE
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  method <- args$method
  run_idx <- as.integer(args$run_idx)
  commit <- args$commit
  date <- args$date
  mode <- args$mode
  threshold <- args$threshold

  random_seed <- args$random_seed
  force_commits <- args$force_commits
  
  output_dir <- args$output_dir
  method_dir <- args$method_dir
}

if (! exists("run_idx", mode="numeric") & mode != "prod") {
  cli_abort(
    "no run index specified",
    x="You must specify which index of the run config to use to set up the model."
  )
} else if (! exists("run_idx", mode="numeric") & mode == "prod") {
  run_idx <- 1
}

if (! exists("random_seed")) {
  random_seed <- NULL
}

if (! exists("threshold")) {
  threshold <- NULL
}

if (! exists("method", mode="character") & ! exists("commit", mode="character") & ! exists("date", mode="character")) {
  cli_abort(c(
    "no method, commit hash, or creation date provided",
    "x"="You must specify which run config to load using at least the {.var method} name."
  ))
}

if (! exists("output_dir")) {
  cli_abort(c(
    "no path to save results provided",
    "x"="You must provide the save path as the variable {.var output_dir}."
  ))
}

if (! mode %in% c("eval", "tune", "final", "prod")) {
  cli_abort(c(
    "unrecognised value {.val {mode}} for {.var mode}",
    x="{.var mode} must be one of {.val {c('eval', 'tune', 'final', 'prod')}}"
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

all_commits <- system2("git", args=c("rev-list", "--max-count=1000", "--abbrev-commit", "--skip=#", "HEAD"),
                       stdout=TRUE)
available_commits <- 
    output_dir |>
    file.path(method) |>
    list.dirs(full.names=FALSE, recursive=FALSE)

if (! exists("commit", mode="character")) {
    latest_idx <- min(which(all_commits %in% available_commits))
    commit <- all_commits[latest_idx]
}

if (! exists("date", mode="character")) {
  all_dates <- list.dirs(file.path(output_dir, method, commit), full.names=FALSE, recursive=FALSE)
  date <- max(all_dates)
}

warn_uncommitted(.stop=force_commits)

output_dir <- file.path(output_dir, method, commit, date)

cli_alert_info("Training a model with line {.strong {run_idx}} of {.strong {mode} scheme} for model {.var {method}-{commit}-{date}}")

## load files ----
metadata <- jsonlite::read_json(file.path(output_dir, "run-metadata.json"))
if (mode != "prod") {
  cv_scheme <- jsonlite::read_json(file.path(output_dir, paste(mode, "scheme.json", sep="-")), simplifyVectors=TRUE)
}

features <- read_rds(file.path(output_dir, "features-labelled.rds"))
unlabelled <- read_rds(file.path(output_dir, "features-unlabelled.rds"))

## data budget ----
if (mode == "prod") {
  train <- features
  test <- unlabelled
} else {
  cv_split <- cv_scheme[[run_idx]]
  train <- filter(features, plant_name_id %in% unlist(cv_split$train))
  test <- filter(features, plant_name_id %in% unlist(cv_split$test))
}

## hyper parameters ----
if (metadata$state == "finalised") {
  hparams <- metadata$hparams
} else {
  hparams <- cv_split$hparams
}

threshold <- hparams$threshold
hparams$threshold <- NULL

## model setup ----
model <- specify_model(target=metadata$target)
model <- finalize_model(model, hparams)

## preprocess data ----
preproc <- specify_recipe(features)
preproc <- prep(preproc, train)

train_proc <- bake(preproc, new_data=NULL)
test_proc <- bake(preproc, new_data=test)

predictor_names <- filter(preproc$var_info, role %in% c("outcome", "predictor"))
train_input <- select(train_proc, all_of(predictor_names$variable))
test_input <- select(test_proc, all_of(predictor_names$variable))

## fit model ----
model_fit <- fit(model, obs ~., data=train_input)

## predict on test set ----
test_preds <- augment(model_fit, new_data=test_input)
if (is.null(threshold) & metadata$target == "threat_status" & mode %in% c("tune", "final")) {
  threshold <- choose_threshold(test_preds$.pred_threatened, test_preds$obs)
  threshold <- threshold$best
} else if (is.null(threshold) & metadata$target == "threat_status") {
  threshold <- 0.5
}

if (metadata$target == "threat_status") {
  test_preds$.pred_class <- ifelse(test_preds$.pred_threatened > threshold, "threatened", "not threatened")

  test_preds$.pred_class <- factor(test_preds$.pred_class, levels=levels(test_preds$obs))
}

## evaluate performance ----
if (mode %in% c("tune", "final")) {
  metrics <- metric_set(roc_auc, mn_log_loss)
} else {
  metrics <- metric_set(accuracy, sensitivity, specificity, j_index)
}

if (metadata$target == "category" & mode != "prod") {
  performance <- metrics(test_preds, truth=obs, estimate=.pred_class, .pred_LC:.pred_CR)
} else if (mode != "prod") {
  performance <- metrics(test_preds, truth=obs, estimate=.pred_class, .pred_threatened, event_level="second")
}

if (metadata$target == "threat_status" & mode %in% c("tune", "final")) {
  performance$threshold <- threshold
}

## generate outputs ----
dir.create(file.path(output_dir, "results"), showWarnings=FALSE)
if (mode == "prod") {
  dir.create(file.path(output_dir, "model"), showWarnings=FALSE)
}

### performance ----
if (mode != "prod") {
  write_csv(performance, file.path(output_dir, "results", paste0(mode, "_idx-", run_idx, "_performance.csv")))
}

### test set predictions ----
if (mode == "eval") {
  write_csv(test_preds, file.path(output_dir, "results", paste0(mode, "_idx-", run_idx, "_test-preds.csv")))
} else if (mode == "prod") {
  write_csv(test_preds, file.path(output_dir, "results", paste0("final-predictions.csv")))
}

if (mode == "prod" & method == "bart") {
  ev_samples <- 
    extract_samples(model_fit, newdata=train_input, type="ev", ids=train$plant_name_id) |>
    mutate(set="labelled") |>
    bind_rows(
      extract_samples(model_fit, newdata=test_input, type="ev", ids=test$plant_name_id) |>
      mutate(set="unlabelled") 
    )

  ppd_samples <- 
    extract_samples(model_fit, newdata=train_input, type="ppd", ids=train$plant_name_id) |>
    mutate(set="labelled") |>
    bind_rows(
      extract_samples(model_fit, newdata=test_input, type="ppd", ids=test$plant_name_id) |>
      mutate(set="unlabelled") 
    )
  
  write_csv(ppd_samples, file.path(output_dir, paste0(name, "-ppd-samples.csv")))
  write_csv(ev_samples, file.path(output_dir, paste0(name, "-ev-samples.csv")))
}

### feature importance ----
if (mode == "eval") {
  importance <- permutation_importance_(model_fit, test_input, parallel=TRUE, .n=5)
  write_csv(importance, file.path(output_dir, "results", paste0(mode, "_idx-", run_idx, "_importance.csv")))
}

### save model ----
if (method == "bart" & mode == "prod") {
  # have to 'touch' the trees to be able to save them for predictions
  invisible(model_fit$fit$fit$state)
}

if (mode == "prod") {
  write_rds(model_fit, file.path(output_dir, "model", paste0(paste(method, commit, date, sep="-"), ".rds")))
  write_rds(preproc, file.path(output_dir, "model", paste0(paste("preprocessor", method, commit, date, sep="-"), ".rds")))
}
