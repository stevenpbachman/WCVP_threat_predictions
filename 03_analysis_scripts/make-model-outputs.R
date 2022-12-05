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

# load predictions for family size
predictions <- read_csv(file.path(model_dir, paste0(name, "-final-predictions.csv")),
                        show_col_types=FALSE)

family_counts <- count(predictions, family)
show_families <- family_counts[family_counts$n > family_size, ]$family

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
  filter(family %in% show_families) |>
  plot_performance_bars("family", "cv", metric="accuracy") +
  facet_wrap(~cv)

climate_perf_plot <- plot_performance_bars(climate_perf, "climate_description", "cv",
                                           metric="j_index")
lifeform_perf_plot <- plot_performance_bars(lifeform_perf, "humphreys_lifeform", "cv",
                                            metric="j_index")

## join into a single plot ----
(perf_plot <- 
  ((overall_perf_plot / (lifeform_perf_plot + remove_xaxis()) / 
     (climate_perf_plot + labs(x="J index")) + plot_layout(heights=c(4, 4, 6))) | 
     (family_perf_plot + labs(x="Accuracy"))) +
 plot_layout(guides="collect") & 
  theme(legend.position="bottom") &
  guides(color="none", fill=guide_legend(title="")))

ggsave(file.path(output_dir, paste0(name, "-performance-plot.png")), perf_plot,
       height=20, width=10)

## region-wise performance map ----
distributions <-
  rWCVPdata::wcvp_distributions |>
  filter(extinct + introduced + location_doubtful == 0) |>
  select(plant_name_id, area_code_l3)

country_perf <-
  distributions |>
  left_join(test_preds, by="plant_name_id") |>
  disaggregate_performance(cv, area_code_l3, metrics=metrics)

(tss_map <-
  country_perf |>
  filter(.metric == "j_index",
         cv == "Random CV") |>
  rename(TSS=.estimate) |>
  plot_map(TSS, .proj="moll", .points=TRUE) +
    scico::scale_fill_scico(
      palette="cork",
      na.value="grey80",
      name="TSS",
      limits=c(-1, 1)
    ) +
    scico::scale_colour_scico(
      palette="cork",
      na.value="grey80",
      name="TSS",
      limits=c(-1, 1)
    ))

(sens_map <-
    country_perf |>
    filter(.metric == "sensitivity",
           cv == "Random CV") |>
    rename(Sensitivity=.estimate) |>
    plot_map(Sensitivity, .proj="moll", .points=TRUE) +
    scico::scale_fill_scico(
      palette="oslo",
      na.value="grey80",
      name="Sensitivity",
      direction=-1,
      limits=c(0, 1)
    ) +
    scico::scale_colour_scico(
      palette="oslo",
      na.value="grey80",
      name="Sensitivity",
      direction=-1,
      limits=c(0, 1)
    ))

(spec_map <-
    country_perf |>
    filter(.metric == "specificity",
           cv == "Random CV") |>
    rename(Specificity=.estimate) |>
    plot_map(Specificity, .proj="moll", .points=TRUE))

ggsave(file.path(output_dir, paste0(name, "-performance-map.png")), perf_map)

test_preds |>
  filter(cv == "Random CV", obs == "threatened") |>
  inner_join(distributions, by="plant_name_id") |>
  group_by(area_code_l3) |>
  summarise(avg_size=log(median(L3_count, na.rm=T))) |>
  plot_map(avg_size, .proj="moll", .points=TRUE)

## richness map ----
distributions |>
  group_by(area_code_l3) |>
  summarise(`Species richness`=n_distinct(plant_name_id)) |>
  plot_map(`Species richness`, .proj="moll", .points=TRUE) +
  scico::scale_fill_scico(
    palette="bamako",
    na.value="grey80",
    name="Species richness",
    direction=-1
  ) +
  scico::scale_colour_scico(
    palette="bamako",
    na.value="grey80",
    name="Species richness",
    direction=-1
  )

## assessments map ----
distributions |>
  filter(plant_name_id %in% predictions[!is.na(predictions$obs),]$plant_name_id) |>
  group_by(area_code_l3) |>
  summarise(`Species richness`=n_distinct(plant_name_id)) |>
  plot_map(`Species richness`, .proj="moll", .points=TRUE) +
  scico::scale_fill_scico(
    palette="oslo",
    na.value="grey80",
    name="Assessed species",
    direction=-1
  ) +
  scico::scale_colour_scico(
    palette="oslo",
    na.value="grey80",
    name="Assessed species",
    direction=-1
  )

# interpretation ----

# variable importance
importance <- read_csv(file.path(model_dir, paste0(name, "-random-importance.csv")),
                       show_col_types=FALSE)

(importance_plot <- 
    importance |>
    group_by(variable) |>
    median_qi(.width=c(0.69, 0.87, 0.95)) |>
    mutate(variable=reorder(variable, mean_decrease_accuracy)) |>
    ggplot(mapping=aes(x=mean_decrease_accuracy, y=variable, xmin=.lower, xmax=.upper)) +
    geom_pointinterval(interval_colour="grey50") +
    geom_vline(xintercept=0, linetype=2, size=1, color="red") +
    labs(x="Mean decrease in accuracy", y=""))

ggsave(file.path(output_dir, paste0(name, "-importance-plot.png")), importance_plot,
       height=10, width=6)

# predictions ----
# load posterior samples if they're there
ppd_file <- list.files(model_dir, pattern=paste0(name, "-ppd-samples"))
if (length(ppd_file) > 0) {
  pred_data <- predictions
  
  predictions <- 
    read_csv(file.path(model_dir, paste0(name, "-ppd-samples.csv")),
             show_col_types=FALSE) |>
    pivot_longer(cols=c(-set, -plant_name_id), names_to=".draw", values_to="threatened") |>
    mutate(.draw=as.integer(str_remove(.draw, "V"))) |>
    bind_cols(
      read_csv(file.path(model_dir, paste0(name, "-ev-samples.csv")),
               show_col_types=FALSE) |>
      pivot_longer(cols=c(-set, -plant_name_id), names_to=".draw", values_to="p_threatened") |>
      select(p_threatened)
    )
  
  predictions <-
    predictions |>
    left_join(
      pred_data, 
      by=c("plant_name_id", "set")
    )
  
  rm(list=c("pred_data"))
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
  filter(family %in% show_families) |>
  plot_threat(family)

(pred_plots <-
  ((overall_preds / lifeform_preds / climate_preds) + 
     plot_layout(heights=c(1, 4, 8)) | family_preds) +
  plot_layout(guides="collect",
              heights=c(1, 1)) &
    theme(legend.position="bottom"))

ggsave(file.path(output_dir, paste0(name, "-threat-plots.png")), pred_plots, 
       height=20, width=10)

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
  plot_map(threatened, .proj="moll", .points=TRUE) +
  scico::scale_fill_scico(
    palette="lajolla",
    na.value="grey80",
    limits=c(0, 1),
    labels=scales::label_percent(),
    name="Predicted threatened"
  ) +
  scico::scale_colour_scico(
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

