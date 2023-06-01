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

#method <- "bart"
#output_dir <- "output/review1/latest_bart_rl2022_2_all/"
#model_dir <- "output/review1/latest_bart_rl2022_2_all/"
#method_dir <- "04_model_definitions"
#predictor_file <- "output/review1/latest_bart_rl2022_2_all/predictors-angiosperm-20230525-155316.csv"
#random_seed <- 1989
#target="threat_status" # or category

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(tidyverse)   # packages for data handling
shhlibrary(tidymodels)  # packages for model training and evaluation
shhlibrary(cli)         # nice formatting for CLI
shhlibrary(git2r)

source("R/utility-functions.R")

# CLI ----
cli_h1("Setting up an extinction prediction model training run")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="output/review1/latest_bart_rl2022_2_all/",
    method_dir="04_model_definitions",
    target="threat_status", # or category
    random_seed=1989,
    force_commits=TRUE
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  method <- args$method
  predictor_file <- args$predictor_file
  target <- args$target
  
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

if (! target %in% c("threat_status", "category")) {
  cli_abort(c(
    "unrecognised target",
    x="Possible targets are the full IUCN RL category ({.val category}) or the threat status ({.val threat_status})"
  ))
}

source(file.path(method_dir, paste0(method, ".R")))

warn_uncommitted(.stop=force_commits)

latest_hash <- system2("git", args=c("rev-list", "--max-count=1", "--abbrev-commit", "--skip=#", "HEAD"),
                       stdout=TRUE)
now <- format(Sys.time(), "%Y%m%d-%H%M%S")

output_dir <- file.path(output_dir, method, now) #latest_hash
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

cli_alert_info("Setting up a {.strong {method}} model using {.file {predictor_file}}")
cli_alert_info("Saving run config to {.file {output_dir}}")

# set up run metadata ----
run_meta <- list(
  commit=paste0("https://github.com/stevenpbachman/WCVP_threat_predictions/commit/", latest_hash),
  date=Sys.time(),
  state="initialised",
  model=method,
  target=target,
  seed=random_seed,
  predictors=predictor_file
)

write(jsonlite::toJSON(run_meta, auto_unbox=TRUE), file.path(output_dir, "run-metadata.json"))

cli_alert_info("predictor file loading next...")
# load predictors ----
predictors <- readr::read_csv(predictor_file, show_col_types=FALSE, progress=FALSE)


cli_alert_info("predictor file loaded")

# recode categories
rl_codes_map <- c("Critically Endangered" = "CR",
                 "Data Deficient" = "DD",
                 "Endangered" = "EN",
                 "Extinct" = "EX",
                 "Extinct in the Wild" = "EW",
                 "Least Concern" = "LC",
                 "Lower Risk/conservation dependent" = "LR/cd",
                 "Lower Risk/least concern" = "LR/lc",
                 "Lower Risk/near threatened" = "LR/nt",
                 "Near Threatened" = "NT",
                 "Vulnerable" = "VU"                       
                 )

predictors$category <- recode(predictors$category, !!! rl_codes_map)

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
} else if (target == "category") {
  labelled$obs <- factor(labelled$category, levels=c("LC", "NT", "VU", "EN", "CR"))
}

unlabelled$obs <- factor(NA, levels=levels(labelled$obs))
write_rds(labelled, file.path(output_dir, "features-labelled.rds"))
write_rds(unlabelled, file.path(output_dir, "features-unlabelled.rds"))

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

final_cv <- vfold_cv(labelled, v=5)

# save evaluation cv scheme ----
random_eval <- insert_items(get_cv_ids(random_cv), list("type"="random"))
family_eval <- insert_items(get_cv_ids(family_cv), list("type"="family"))
eval_scheme <- c(random_eval, family_eval)
write(jsonlite::toJSON(as.list(eval_scheme), auto_unbox=TRUE), file.path(output_dir, "eval-scheme.json"))

# save tune cv scheme ----
# need to expand tuning scheme so each hyperparameter is evaluated on each fold
if (exists("hparam_grid")) {
  random_tune_ids <- map2(random_cv$inner_resamples, random_cv$id, ~insert_items(get_cv_ids(.x), list("outer"=.y)))
  random_tune_ids <- do.call(c, random_tune_ids)
  random_tune_ids <- insert_items(random_tune_ids, list("type"="random"))
  
  family_tune_ids <- map2(family_cv$inner_resamples, family_cv$id, ~insert_items(get_cv_ids(.x), list("outer"=.y)))
  family_tune_ids <- do.call(c, family_tune_ids)
  family_tune_ids <- insert_items(family_tune_ids, list("type"="family"))

  hparams <- jsonlite::toJSON(hparam_grid)
  hparam_list <- jsonlite::fromJSON(hparams, simplifyDataFrame=FALSE)
  hparam_random_cv <- do.call(c, lapply(hparam_list, function(x) insert_items(random_tune_ids, list("hparams"=x))))
  hparam_family_cv <- do.call(c, lapply(hparam_list, function(x) insert_items(family_tune_ids, list("hparams"=x))))
  
  tune_scheme <- c(hparam_random_cv, hparam_family_cv)
  write(jsonlite::toJSON(tune_scheme, auto_unbox=TRUE), file.path(output_dir, "tune-scheme.json"))
}

# save final cv scheme ----
# only need this if there are hyperparameters to tune
if (exists("hparam_grid")) {
  final_tune_ids <- insert_items(get_cv_ids(final_cv), list("type"="random"))

  hparams <- jsonlite::toJSON(hparam_grid)
  hparam_list <- jsonlite::fromJSON(hparams, simplifyDataFrame=FALSE)
  final_tune_scheme <- do.call(c, lapply(hparam_list, function(x) insert_items(final_tune_ids, list("hparams"=x))))
  
  write(jsonlite::toJSON(final_tune_scheme, auto_unbox=TRUE), file.path(output_dir, "final-scheme.json"))
}
