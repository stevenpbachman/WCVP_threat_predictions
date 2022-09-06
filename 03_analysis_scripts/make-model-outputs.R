#' Script to train, evaluate, and make predictions of extinction risk

#' EXPECTED INPUTS:
#'  - `model_name`: the name of the model to generate outputs for, matches partially
#'                  so `bart` will make outputs for the latest bart model, whereas
#'                  `bart-8b27a1e209ec72aee44e440e87a97c58444046ef-20220830-203114`
#'                  will make outputs for that specific bart model.
#'  - `output_dir`: path to a directory to save outputs to.
#'  - `model_dir`: path where the model and related files are saved.
#' 
#' EXAMPLE CLI:
#'  Rscript 03_analysis_scripts/make-model-outputs.R --model_name=bart
#' 
#' EXAMPLE SOURCE:
#'  output_dir <- "output"
#'  model_dir <- "output"
#'  model_name <- "bart"
#'  
#'  source("03_analysis_scripts/make-model-outputs.R")
#' 

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(tidyverse)   # packages for data handling
shhlibrary(lubridate)   # handle date time conversions
shhlibrary(yardstick)   # functions for calculating performance metrics
shhlibrary(cli)         # nice formatting for CLI
shhlibrary(patchwork)   # joining multiple plots together
shhlibrary(ggdist)      # ggplot functions for visualising distributions
shhlibrary(sf)          # handle and plot shape files

source("R/plot-functions.R")
source("R/utility-functions.R")

# CLI ----
cli_h1("Running extinction risk prediction method")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="output",
    model_dir="output"
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  model_name <- args$model_name
  output_dir <- args$output_dir
  model_dir <- args$model_dir
}

if (! exists("model_name", mode="character")) {
  cli_abort(c(
    "no model name provided",
    "x"="You must provide the name of a trained model to make outputs for as {.var model_name}."
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

dir.create(output_dir, showWarnings=FALSE)

models <- list.files(model_dir, pattern=paste0(model_name, ".*.rds"))

name <- str_remove(models[length(models)], ".rds")

name_parts <- str_split(name, "-")[[1]]
model_hash <- name_parts[length(name_parts) - 2]
model_date <- name_parts[length(name_parts) - 1]
model_time <- name_parts[length(name_parts)]

model_datetime <- lubridate::as_datetime(paste(model_date, model_time))

cli_alert_info("Generating outputs for {.strong {model_name}} model trained on {.strong {model_datetime}} using code from {.url {paste0('https://github.com/stevenpbachman/WCVP_threat_predictions/commit/', model_hash)}}")
cli_alert_info("Saving outputs to {.file {output_dir}}")

# performance ----
## load data ----
random_performance <- read_csv(file.path(model_dir, paste0(name, "-random-performance.csv")),
                               show_col_types=FALSE)
family_performance <- read_csv(file.path(model_dir, paste0(name, "-family-performance.csv")),
                               show_col_types=FALSE)

random_test_preds <- read_csv(file.path(model_dir, paste0(name, "-random-test-preds.csv")),
                              show_col_types=FALSE)
family_test_preds <- read_csv(file.path(model_dir, paste0(name, "-family-test-preds.csv")),
                              show_col_types=FALSE)

## overall performance ----
overall_performance <-
  random_performance |>
  mutate(cv="Random CV") |>
  bind_rows(
    family_performance |>
      mutate(cv="Family-wise block CV")
  ) |>
  mutate(cv=factor(cv, levels=c("Random CV", "Family-wise block CV"))) |>
  group_by(cv, .metric) |>
  mean_hdci(.estimate)

(overall_perf_plot <- ggplot(data=overall_performance,
                             mapping=aes(x=.estimate, y=.metric, xmax=.upper, xmin=.lower,
                                         colour=cv)) +
  geom_pointrange(size=1, position=position_dodge(width=0.4)) +
  labs(x="Value", y="") +
  lims(x=c(0, 1)))

## disaggregated performance ----
metrics <- metric_set(accuracy, sensitivity, specificity, j_index)

test_preds <- 
  random_test_preds |>
  mutate(cv="Random CV") |>
  bind_rows(
    family_test_preds |>
      mutate(cv="Family-wise block CV")
  ) |>
  mutate(cv=factor(cv, levels=c("Random CV", "Family-wise block CV"))) |>
  mutate(obs=factor(obs, levels=c("not threatened", "threatened"))) |>
  mutate(.pred_class=factor(.pred_class, levels=levels(obs)))

family_perf <- disaggregate_performance(test_preds, cv, family, metrics=metrics)
lifeform_perf <- disaggregate_performance(test_preds, cv, humphreys_lifeform, metrics=metrics)
climate_perf <- disaggregate_performance(test_preds, cv, climate_description, metrics=metrics)

family_perf_plot <- 
  family_perf |>
  plot_performance_bars("family", "cv", metric="j_index") +
  facet_wrap(~cv)

climate_perf_plot <- plot_performance_bars(climate_perf, "climate_description", "cv",
                                           metric="j_index")
lifeform_perf_plot <- plot_performance_bars(lifeform_perf, "humphreys_lifeform", "cv",
                                            metric="j_index")

## join into a single plot ----
(perf_plot <- 
  ((overall_perf_plot / (lifeform_perf_plot + remove_xaxis()) / 
     (climate_perf_plot + labs(x="J index")) + plot_layout(heights=c(4, 4, 6))) | 
     (family_perf_plot + labs(x="J index"))) +
 plot_layout(guides="collect") & 
  theme(legend.position="bottom") &
  guides(color="none", fill=guide_legend(title="")))

ggsave(file.path(output_dir, paste0(name, "-performance-plot.png")), perf_plot,
       height=10, width=10)

## region-wise performance map ----
distributions <-
  rWCVPdata::wcvp_distributions |>
  filter(extinct + introduced + location_doubtful == 0) |>
  select(plant_name_id, area_code_l3)

country_perf <-
  distributions |>
  left_join(test_preds, by="plant_name_id") |>
  disaggregate_performance(cv, area_code_l3, metrics=metrics)

(perf_map <-
  country_perf |>
  filter(.metric == "j_index",
         cv == "Random CV") |>
  rename(`J index`=.estimate) |>
  plot_map(`J index`))

ggsave(file.path(output_dir, paste0(name, "-performance-map.png")), perf_map)

# predictions ----
## load data ----
predictions <- read_csv(file.path(model_dir, paste0(name, "-final-predictions.csv")),
                        show_col_types=FALSE)

# load posterior samples if they're there
ppd_file <- list.files(model_dir, pattern=paste0(name, "-ppd-samples"))
if (length(ppd_file) > 0) {
  pred_samples <- read_csv(file.path(model_dir, paste0(name, "-ppd-samples.csv")),
                           show_col_types=FALSE)
  pred_probs <- read_csv(file.path(model_dir, paste0(name, "-ev-samples.csv")),
                         show_col_types=FALSE)
  
  pred_samples <- 
    pred_samples |>
    pivot_longer(cols=c(-set, -plant_name_id), names_to=".draw", values_to="threatened") |>
    mutate(.draw=as.integer(str_remove(.draw, "V"))) 
  
  pred_probs <- 
    pred_probs |>
    pivot_longer(cols=c(-set, -plant_name_id), names_to=".draw", values_to="p_threatened") |>
    mutate(.draw=as.integer(str_remove(.draw, "V")))
  
  predictions <-
    pred_samples |>
    left_join(pred_probs, by=c("plant_name_id", "set", ".draw")) |>
    left_join(predictions, by=c("plant_name_id", "set"))
}

## plot predictions of % threatened ----
overall_preds <- 
  predictions |>
  mutate(name="Overall") |>
  plot_threat(name) +
  remove_xaxis()

climate_preds <-
  predictions |>
  plot_threat(climate_description)

lifeform_preds <-
  predictions |>
  plot_threat(humphreys_lifeform) +
  remove_xaxis()

family_preds <-
  predictions |>
  plot_threat(family)

(pred_plots <-
  ((overall_preds / lifeform_preds / climate_preds) + 
     plot_layout(heights=c(1, 4, 8)) | family_preds) +
  plot_layout(guides="collect",
              heights=c(1, 1)) &
    theme(legend.position="bottom"))

ggsave(file.path(output_dir, paste0(name, "-threat-plots.png")), pred_plots, 
       height=10, width=10)

## predictions for each WGSRPD L3 region ----
if (".draw" %in% colnames(predictions)) {
  country_preds <-
    distributions |>
    left_join(
      predictions |> select(plant_name_id, .draw, threatened),
      by="plant_name_id"
    ) |>
    group_by(area_code_l3, .draw)
} else {
  country_preds <-
    distributions |>
    left_join(
      predictions |> select(plant_name_id, .pred_class),
      by="plant_name_id"
    ) |>
    mutate(threatened=ifelse(.pred_class == "threatened", 1, 0)) |>
    group_by(area_code_l3)
}

country_preds <-
  country_preds |>
  summarise(threatened=mean(threatened, na.rm=TRUE))

if (".draw" %in% colnames(country_preds)) {
  country_preds <-
    country_preds |>
    median_qi(threatened, na.rm=TRUE)
}

(threat_map <- 
  country_preds |>
  plot_map(threatened) +
  scico::scale_fill_scico(
    palette="lajolla",
    na.value="grey80",
    limits=c(0, 1),
    labels=scales::label_percent(),
    name="Predicted threatened"
  ))

ggsave(file.path(output_dir, paste0(name, "-threat-map.png")), threat_map)

if (".draw" %in% colnames(predictions)) {
  (post_map <- 
    country_preds |>
    mutate(`Posterior width`=.upper - .lower) |>
    plot_map(`Posterior width`) +
    scico::scale_fill_scico(
      palette="lapaz",
      na.value="grey80",
      limits=c(0, 1),
      direction=-1
    ))
  ggsave(file.path(output_dir, paste0(name, "-uncertainty-map.png")), post_map)
}

# interpretation ----

# variable importance
importance <- read_csv(file.path(model_dir, paste0(name, "-random-importance.csv")),
                       show_col_types=FALSE)

importance_plot <-
  ggplot(data=importance, mapping=aes(x=mean_decrease_accuracy, y=variable)) +
  geom_boxplot() +
  labs(x="Mean decrease in accuracy", y="")

ggsave(file.path(output_dir, paste0(name, "-importance-plot.png")), importance_plot)
