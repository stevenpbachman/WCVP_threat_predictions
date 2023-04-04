#' Set up the hyperparameters, cross-validation folds, and pre-processor for 
#' training a model.


#' EXPECTED INPUTS:
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
shhlibrary(dbarts)      # for a bart model if needed
shhlibrary(tidyverse)   # packages for data handling
shhlibrary(tidymodels)  # packages for model training and evaluation
shhlibrary(ROCR)        # evaluate classification thresholds
shhlibrary(multidplyr)  # run operations on chunks of a data frame in parallel
shhlibrary(doParallel)  # set up parallel processing
shhlibrary(git2r)       # run git commands in R
shhlibrary(cli)         # nice formatting for CLI


# CLI ----
cli_h1("Setting up an extinction prediction model training run")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="05_outputs",
    method_dir="04_model_definitions",
    target="threat_status",
    random_seed=1989,
    force_commits=TRUE
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  method <- args$method
  predictor_file <- args$predictor_file
  
  random_seed <- args$random_seed
  force_commits <- args$force_commits
  
  output_dir <- args$output_dir
  method_dir <- args$method_dir
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

warn_uncommitted(.stop=force_commits)

latest_hash <- revparse_single(".", "HEAD")$sha
now <- format(Sys.time(), "%Y%m%d-%H%M%S")
name <- paste(method, latest_hash, now, sep="-")

cli_alert_info("Setting up a {.strong {method}} model using {.file {predictor_file}}")
cli_alert_info("Saving run config to {.file {output_dir}}")

# load predictors ----
predictors <- read_csv(predictor_file, show_col_types=FALSE, progress=FALSE)

# prepare predictors ----
predictors <- filter(predictors, ! category %in% c("EW", "EX"))

old_codes_map <- c("LR/cd"="NT", "LR/lc"="LC", "LR/nt"="NT")
predictors$category <- recode(predictors$category, !!! old_codes_map)

factor_vars <- c("humphreys_lifeform", "climate_description", "genus", "family", "order", "higher_groups")
predictors <- mutate(predictors, across(all_of(factor_vars), ~factor(.x)))

labelled <- filter(predictors, ! is.na(category), category != "DD")
unlabelled <- filter(predictors, is.na(category) | category == "DD")

if (target == "threat_status") {
  labelled$obs <- ifelse(labelled$category %in% c("LC", "NT"), "not threatened", "threatened")
  labelled$obs <- factor(labelled$obs, levels=c("not threatened", "threatened"))  
} else if (target == "categories") {
  labelled$obs <- factor(labelled$categories, levels=c("LC", "NT", "VU", "EN", "CR"))
}

unlabelled$obs <- factor(NA, levels=levels(labelled$obs))

cli_alert_info("Evaluating the method on {.strong {nrow(labelled)}} examples:")

cat_counts <- count(labelled, obs)
cat_counts$p <- cat_counts$n / sum(cat_counts$n)

for (i in seq_along(cat_counts$obs)) {
  cli_text("{cat_counts$obs[i]}: {cat_counts$n[i]} ({scales::percent_format(accuracy=0.1)(cat_counts$p[i])})")
}

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

# this sets up a cluster for mapping chunks of data frames (model evaluation)
if (parallel) {
  cluster <- new_cluster(ncores)
  cluster_library(cluster, packages=c("dplyr", "dbarts", "tidymodels", "ROCR", "vip"))
} else {
  cluster <- NULL
}

# set up model ----
model_spec <- specify_model()
model_recipe <- specify_recipe(labelled)

wf <- 
  workflow() |>
  add_model(model_spec) |>
  add_recipe(model_recipe)

# tune hyperparameters if needed ----
model_params <- extract_parameter_set_dials(wf)
untuned_params <- length(model_params$name)