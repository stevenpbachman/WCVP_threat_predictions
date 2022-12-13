#' Script to create conformal predictions for the output of a model

#' EXPECTED INPUTS:
#'  - `model_name`: the name of the model to generate outputs for, matches partially
#'                  so `bart` will make outputs for the latest bart model, whereas
#'                  `bart-8b27a1e209ec72aee44e440e87a97c58444046ef-20220830-203114`
#'                  will make outputs for that specific bart model.
#'  - `output_dir`: path to a directory to save outputs to.
#'  - `model_dir`: path where the model and related files are saved.
#' 
#' EXAMPLE CLI:
#'  Rscript 03_analysis_scripts/conformal-predictions.R --model_name=bart
#' 
#' EXAMPLE SOURCE:
#'  output_dir <- "output"
#'  model_dir <- "output"
#'  model_name <- "bart"
#'  
#'  source("03_analysis_scripts/conformal-predictions.R")
#' 

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(tidyverse)   # packages for data handling
shhlibrary(lubridate)   # handle date time conversions
shhlibrary(cli)         # nice formatting for CLI
shhlibrary(patchwork)   # joining multiple plots together
shhlibrary(ggdist)      # ggplot functions for visualising distributions
shhlibrary(sf)          # handle and plot shape files

source("R/plot-functions.R")
source("R/utility-functions.R")
source("R/conformal-functions.R")

# CLI ----
cli_h1("Generating prediction sets")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="output",
    model_dir="output",
    family_size=0
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  model_name <- args$model_name
  output_dir <- args$output_dir
  model_dir <- args$model_dir
  family_size <- args$family_size
}

if (! exists("model_name", mode="character")) {
  cli_abort(c(
    "no model name provided",
    "x"="You must provide the name of a trained model to make prediction sets for as {.var model_name}."
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
    "no path to model results provided",
    "x"="You must provide the path to where model results are saved as the variable {.var model_dir}."
  ))
}

if (! exists("family_size", mode="integer")) {
  family_size <- 0
}

dir.create(output_dir, showWarnings=FALSE)

models <- list.files(model_dir, pattern=paste0(model_name, ".*.rds"))

name <- str_remove(models[length(models)], ".rds")

name_parts <- str_split(name, "-")[[1]]
model_hash <- name_parts[length(name_parts) - 2]
model_date <- name_parts[length(name_parts) - 1]
model_time <- name_parts[length(name_parts)]

model_datetime <- lubridate::as_datetime(paste(model_date, model_time))

cli_alert_info("Generating predictin sets for {.strong {model_name}} model trained on {.strong {model_datetime}} using code from {.url {paste0('https://github.com/stevenpbachman/WCVP_threat_predictions/commit/', model_hash)}}")
cli_alert_info("Saving outputs to {.file {output_dir}}")

# performance ----
## load data ----
test_preds <- read_csv(file.path(model_dir, paste0(name, "-random-test-preds.csv")),
                              show_col_types=FALSE)

# load predictions
predictions <- read_csv(file.path(model_dir, paste0(name, "-final-predictions.csv")),
                        show_col_types=FALSE)


# calibrate ----

scores <- calculate_scores(test_preds, label_col="obs")
thresholds <- 
  scores |>
  group_by(id) |>
  summarise(qhat=calculate_threshold(.score, alpha=0.1))

overall_threshold <- mean(thresholds$qhat)

# check calibration ----

## marginal ----
marginal_coverage <- 
  scores |>
  nest_by(id) |>
  mutate(coverage=list(estimate_coverage(data, alpha=0.1, times=1000))) |>
  select(-data) |>
  unnest(coverage)

write_csv(marginal_coverage, file.path(output_dir, paste0(name, "-marginal-coverage.csv")))

## class conditional ----
conditional_coverage <- 
  scores |>
  nest_by(id) |>
  mutate(coverage=list(estimate_coverage(data, alpha=0.1, times=1000, cls_conditional=TRUE))) |>
  select(-data) |>
  unnest(coverage)

write_csv(conditional_coverage, file.path(output_dir, paste0(name, "-conditional-coverage.csv")))

# generate prediction sets ----
prediction_sets <- make_prediction_set(predictions, overall_threshold)

write_csv(prediction_sets, file.path(output_dir, paste0(name, "-prediction-sets.csv")))
