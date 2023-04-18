#' Update model config files after tuning or evaluation

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
shhlibrary(tidyverse)
shhlibrary(cli)         # nice formatting for CLI

source("R/utility-functions.R")
# CLI ----
cli_h1("Setting up an extinction prediction model training run")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="05_outputs"
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  method <- args$method
  commit <- args$commit
  date <- args$date
  
  output_dir <- args$output_dir
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

available_methods <- list.dirs(output_dir, full.names=FALSE, recursive=FALSE)
if (! method %in% available_methods) {
    cli_abort(c(
        "no config file found for method {.var {method}}",
        x="Config files available for {.var {available_methods}}. Are you sure you ran {.var setup-model.R}?"
    ))
}

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

output_dir <- file.path(output_dir, method, commit, date)

cli_alert_info("Updating config for model {.var {method}-{commit}-{date}}")

## load metadata ----
metadata <- jsonlite::read_json(file.path(output_dir, "run-metadata.json"))

eval_scheme <- jsonlite::read_json(file.path(output_dir, "eval-scheme.json"))
final_scheme <- jsonlite::read_json(file.path(output_dir, "final-scheme.json"))

updated_metadata <- metadata

## update metadata ----
if (metadata$state == "initialised" & file.exists(file.path(output_dir, "tune-scheme.json"))) {
  ### after tuning ----
  tune_scheme <- jsonlite::read_json(file.path(output_dir, "tune-scheme.json"))
  tune_files <- list.files(file.path(output_dir, "results"), pattern="tune_.*_performance.csv", full.names=TRUE)
  if (length(tune_files) != length(tune_scheme)) {
    cli_abort(
      "number of tuning results different from number of tuning folds",
      x="it looks like some tuning results files are missing, have you finished tuning the model?"
    )
  }

  tune_results <- 
    tune_files |>
    read_csv(id="fname") |>
    mutate(run_id=as.integer(str_extract(fname, "\\d+(?=_performance)"))) |>
    select(-fname)

  tune_info <- bind_cols(
    map_dfr(tune_scheme, ~as_tibble(.x[c("type", "outer", "name")])),
    map_dfr(tune_scheme, ~as_tibble(.x$hparams))
  )

  hparam_names <- names(tune_scheme[[1]]$hparams)
  tune_info$run_id <- seq(from=1, to=length(tune_scheme))
  
  tune_results <- 
    tune_results |>
    left_join(tune_info, by="run_id", relationship="many-to-one")

  summary_cols <- c("type", "outer", ".metric", ".estimator", hparam_names)
  tune_summary <-
      tune_results |>
      group_by(across(all_of(summary_cols))) |>
      summarise(
        mean=mean(.estimate, na.rm=TRUE),
        n=n(),
        std_err=sd(.estimate, na.rm=TRUE) / sqrt(n),
        .groups="drop"
      )

  best_hparams <- tune_summary |>
    filter(.metric == "roc_auc") |>
    group_by(type, outer) |>
    slice_max(mean, n=1) |>
    ungroup() |>
    select(-.metric, -.estimator, -n, -mean, -std_err) |>
    nest_by(type, outer, .key="hparams") |>
    mutate(hparams=list(as.list(hparams)))

  updated_eval_scheme <-
    eval_scheme |>
    map_dfr(~tibble(name=.x$name, type=.x$type, train=list(.x$train), test=list(.x$test))) |>
    left_join(best_hparams, relationship="one-to-one", by=c("name"="outer", "type")) |>
    jsonlite::toJSON(auto_unbox=TRUE)

  updated_metadata$state <- "tuned"

  write(updated_eval_scheme, file.path(output_dir, "eval-scheme.json"))
  write(jsonlite::toJSON(updated_metadata, auto_unbox=TRUE), file.path(output_dir, "run-metadata.json"))
  write_csv(tune_summary, file.path(output_dir, "hparam-tuning-summary.csv"))
} else if (metadata$state %in% c("initialised", "tuned")) {
  ### after evaluations ----
  eval_files <- list.files(file.path(output_dir, "results"), pattern="eval_.*_performance.csv", full.names=TRUE)
  if (length(eval_files) != length(eval_scheme)) {
    cli_abort(
      "number of evaluation results different from number of evaluation folds",
      x="it looks like some evaluation results files are missing, have you finished evaluating the model?"
    )
  }

  eval_results <- 
    eval_files |>
    read_csv(id="fname") |>
    mutate(run_id=as.integer(str_extract(fname, "\\d+(?=_performance)"))) |>
    select(-fname)

  eval_info <- bind_cols(
    map_dfr(eval_scheme, ~as_tibble(.x[c("type", "name")])),
    map_dfr(eval_scheme, ~as_tibble(.x$hparams))
  )
  eval_info$run_id <- seq(from=1, to=length(eval_scheme))
  
  eval_results <- 
    eval_results |>
    left_join(eval_info, by="run_id", relationship="many-to-one")

  updated_metadata$state <- "evaluated"

  write(jsonlite::toJSON(updated_metadata, auto_unbox=TRUE), file.path(output_dir, "run-metadata.json"))
  write_csv(eval_results, file.path(output_dir, "eval-performance-summary.csv"))
} else if (metadata$state == "evaluated" & file.exists(file.path(output_dir, "final-scheme.json"))) {
  ### after final tuning ----
  final_scheme <- jsonlite::read_json(file.path(output_dir, "final-scheme.json"))
  final_files <- list.files(file.path(output_dir, "results"), pattern="final_.*_performance.csv", full.names=TRUE)
  if (length(final_files) != length(final_scheme)) {
    cli_abort(
      "number of final tuning results different from number of final tuning folds",
      x="it looks like some tuning results files are missing, have you finished tuning the final model?"
    )
  }

  final_results <- 
    final_files |>
    read_csv(id="fname") |>
    mutate(run_id=as.integer(str_extract(fname, "\\d+(?=_performance)"))) |>
    select(-fname)

  final_info <- bind_cols(
    map_dfr(final_scheme, ~as_tibble(.x[c("type", "name")])),
    map_dfr(final_scheme, ~as_tibble(.x$hparams))
  )
  final_info$run_id <- seq(from=1, to=length(final_scheme))

  hparam_names <- names(tune_scheme[[1]]$hparams)
  
  final_results <- 
    final_results |>
    left_join(final_info, by="run_id", relationship="many-to-one")

  summary_cols <- c("type", "outer", ".metric", ".estimator", hparam_names)
  final_summary <-
    final_results |>
    group_by(across(all_of(hparam_names))) |>
    summarise(
      mean=mean(.estimate, na.rm=TRUE),
      n=n(),
      std_err=sd(.estimate, na.rm=TRUE) / sqrt(n),
      .groups="drop"
    )

  final_hparams <- 
    final_summary |>
    filter(.metric == "roc_auc") |>
    group_by(type) |>
    slice_max(mean, n=1) |>
    ungroup() |>
    select(-type, -.metric, -.estimator, -n, -mean, -std_err) |>
    as.list()

  updated_metadata$hparams <- final_hparams
  updated_metadata$state <- "finalised"

  write(jsonlite::toJSON(updated_metadata, auto_unbox=TRUE), file.path(output_dir, "run-metadata.json"))
  write_csv(final_summary, file.path(output_dir, "final-tuning-summary.csv"))
} else if (metadata$state == "finalised") {
  updated_metadata$state <- "done"
  updated_metadata$final_model <- paste0(paste(method, commit, date, sep="-"), ".rds")
  write(jsonlite::toJSON(updated_metadata, auto_unbox=TRUE), file.path(output_dir, "run-metadata.json"))
} else {
  cli_abort(
    "unrecognised modelling state",
    x="the state of the model, as recorded in the {.file run-metadata.json}, must be one of 'initialised', 'tuned', 'evaluated', or 'finalised'."
  )
}
