library(tidyverse)
library(lubridate)
library(yardstick)
library(cli)
library(patchwork)
library(ggdist)
library(sf)
library(broom)
library(writexl)
library(multidplyr)

source("R/plot-functions.R")
source("R/utility-functions.R")

# setup ----
library(extrafont)
font_import(prompt=FALSE)
loadfonts(device="all")
theme_set(theme_gray(base_family="sans"))
family_size <- 3500

home_dir <- Sys.getenv("HOME")
model_dir <- file.path(home_dir, "projects/rbgk/users_area/bwalker/WCVP_threat_predictions")
name <- "bart-948d0564e4f598b192b326e49a3ddaf5646ddc5a-20221126-221001"
print(model_dir)

# load data ----
random_performance <- read_csv(file.path(model_dir, paste0(name, "-random-performance.csv")),
                               show_col_types=FALSE)
family_performance <- read_csv(file.path(model_dir, paste0(name, "-family-performance.csv")),
                               show_col_types=FALSE)

random_test_preds <- read_csv(file.path(model_dir, paste0(name, "-random-test-preds.csv")),
                              show_col_types=FALSE)
family_test_preds <- read_csv(file.path(model_dir, paste0(name, "-family-test-preds.csv")),
                              show_col_types=FALSE)

pred_data <- read_csv(file.path(model_dir, paste0(name, "-final-predictions.csv")), show_col_type=FALSE)

predictions <- 
  read_csv(file.path(model_dir, paste0(name, "-ppd-samples.csv")),
           show_col_types=FALSE) |>
  pivot_longer(cols=c(-set, -plant_name_id), names_to=".draw", 
               values_to="threatened", names_prefix="V", 
               names_transform=as.integer)

probabilities <- 
  read_csv(file.path(model_dir, paste0(name, "-ev-samples.csv")),
           show_col_types=FALSE) |>
  pivot_longer(cols=c(-set, -plant_name_id), names_to=".draw", 
               values_to="threatened", names_prefix="V", 
               names_transform=as.integer)
  

importance <- read_csv(file.path(model_dir, paste0(name, "-random-importance.csv")),
                       show_col_types=FALSE)
# performance ----

## overall ----
overall_performance <-
  random_performance |>
  mutate(cv="Random CV") |>
  bind_rows(
    family_performance |>
      mutate(cv="Family-wise block CV")
  ) |>
  mutate(cv=factor(cv, levels=c("Random CV", "Family-wise block CV"))) |>
  mutate(.metric=str_to_title(.metric),
         .metric=recode(.metric, J_index="TSS")) |>
  group_by(cv, .metric) |>
  mean_hdci(.estimate)

fmt <- scales::label_number(accuracy=0.001)
overall_performance |>
  mutate(blah=glue::glue("{fmt(.estimate)}\n[{fmt(.lower)}, {fmt(.upper)}]")) |>
  select(cv, .metric, blah) |>
  pivot_wider(names_from=.metric, values_from=blah) |>
  knitr::kable()

## disaggregated ----
metrics <- metric_set(accuracy, sensitivity, specificity, j_index)

family_counts <- count(pred_data, family)
show_families <- family_counts[family_counts$n > family_size, ]$family

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

wgsrpd3 <- rWCVPdata::wgsrpd3
classify_zone <- function(x) {
    centroid <- st_centroid(x)
    centroid_coords <- st_coordinates(centroid)
    lat <- centroid_coords[,2]
    case_when(
      abs(lat) <= 23.3 ~ "Tropical",
      abs(lat) > 23.3 & abs(lat) <= 30 ~ "Subtropical",
      abs(lat) > 30 & abs(lat) <= 66.3 ~ "Temperate",
      TRUE ~ "Polar"
    )
}

zones <- 
  rWCVPdata::wgsrpd3 |>
  mutate(zone=classify_zone(geometry)) |>
  st_drop_geometry() |>
  select(LEVEL3_COD, zone)

species_zones <-
  rWCVPdata::wcvp_distributions |>
  filter(plant_name_id %in% pred_data$plant_name_id) |>
  filter(location_doubtful + extinct + introduced == 0) |>
  left_join(zones, by=c("area_code_l3"="LEVEL3_COD"))

zone_perf <- test_preds |>
    select(cv, plant_name_id, obs, .pred_class) |>
    left_join(species_zones, by="plant_name_id") |>
    disaggregate_performance(cv, zone, metrics=metrics)

## plots ----
cv_colours <- c(
  "Random CV"="#024b7a",
  "Family-wise block CV"="#44b7c2"
)

family_perf_plot <- 
  family_perf |>
  filter(family %in% show_families) |>
  plot_performance_bars("family", "cv", metric="accuracy") +
  scale_colour_manual(values=cv_colours)

zone_perf_plot <- 
  zone_perf |>
  plot_performance_bars("zone", "cv", metric="j_index") +
  scale_colour_manual(values=cv_colours)

lifeform_perf_plot <- 
  lifeform_perf |>
  replace_na(list(humphreys_lifeform="unknown")) |>
  mutate(humphreys_lifeform=str_to_sentence(humphreys_lifeform)) |>
  plot_performance_bars("humphreys_lifeform", "cv", metric="j_index") +
  scale_colour_manual(values=cv_colours)

climate_perf_plot <- 
  climate_perf |>
  replace_na(list(climate_description="Unknown")) |>
  plot_performance_bars("climate_description", "cv", metric="j_index") +
  scale_colour_manual(values=cv_colours)
  

perf_plots_climate <- (((lifeform_perf_plot + labs(x="")) / 
     (climate_perf_plot + labs(x="TSS")) +
    plot_layout(heights=c(5, 8))) | 
     (family_perf_plot + labs(x="Accuracy"))) +
 plot_layout(guides="collect", heights=c(1, 1), ) +
  plot_annotation(tag_levels="A", ) & 
  theme(legend.position="bottom", 
        panel.background=element_blank(),
        panel.grid.major.x=element_line(colour="grey80", linetype=2),
        panel.grid.major.y=element_blank(),
        text=element_text(family="DejaVu Sans")) &
  guides(fill="none", colour=guide_legend(title=""))
ggsave("output/bart-performance-disaggregated-climate.svg", perf_plots_climate, height=10, width=7)

perf_plots_zone <- (((lifeform_perf_plot + labs(x="")) / 
     (zone_perf_plot + labs(x="TSS")) +
    plot_layout(heights=c(5, 4))) | 
     (family_perf_plot + labs(x="Accuracy"))) +
 plot_layout(guides="collect", heights=c(1, 1), ) +
  plot_annotation(tag_levels="A", ) & 
  theme(legend.position="bottom", 
        panel.background=element_blank(),
        panel.grid.major.x=element_line(colour="grey80", linetype=2),
        panel.grid.major.y=element_blank(),
        text=element_text(family="DejaVu Sans")) &
  guides(fill="none", colour=guide_legend(title=""))
ggsave("output/bart-performance-disaggregated-zones.svg", perf_plots_zone, height=10, width=7)

## maps ----
distributions <-
  rWCVPdata::wcvp_distributions |>
  filter(extinct + introduced + location_doubtful == 0) |>
  select(plant_name_id, area_code_l3)

country_perf <-
  distributions |>
  left_join(test_preds, by="plant_name_id") |>
  disaggregate_performance(cv, area_code_l3, metrics=metrics)

tss_map <-
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
    )

sens_map <-
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
    )

spec_map <-
    country_perf |>
    filter(.metric == "specificity",
           cv == "Random CV") |>
    rename(Specificity=.estimate) |>
    plot_map(Specificity, .proj="moll", .points=TRUE) +
    scico::scale_fill_scico(
      palette="acton",
      na.value="grey80",
      name="Specificity",
      direction=-1,
      limits=c(0, 1)
    ) +
    scico::scale_colour_scico(
      palette="acton",
      na.value="grey80",
      name="Specificity",
      direction=-1,
      limits=c(0, 1)
    )

vertical <- tss_map / sens_map / spec_map + plot_annotation(tag_levels="A") & theme(text=element_text(family="DejaVu Sans"))
ggsave("output/bart-performance-maps-vertical.svg", vertical,
       height=9, width=6)

# feature importance ----

full_importance_plot <-
  importance |>
    group_by(variable) |>
    median_qi(.width=0.95) |>
    mutate(variable=reorder(variable, mean_decrease_accuracy)) |>
    ggplot(mapping=aes(x=mean_decrease_accuracy, y=variable, xmin=.lower, xmax=.upper, colour=mean_decrease_accuracy)) +
    geom_pointinterval(colour="grey40", show.legend=FALSE) +
    geom_vline(xintercept=0, linetype=2, size=1, color="red") +
    labs(x="Mean decrease in accuracy", y="") +
  theme(panel.background=element_blank(),
        panel.grid.major.x=element_line(linetype=2, colour="grey80"),
        panel.grid.major.y=element_blank(),
        text=element_text(family="DejaVu Sans"))

ggsave("output/bart-feature-importance-full.svg", full_importance_plot, height=10, width=6)

grouped_importance_plot <-
importance |>
  mutate(fold=c(rep(1, 88), rep(2, 88), 
                rep(3, 88), rep(4, 88), rep(5, 88))) |>
  filter(! variable %in% c("genus", "family", "order", "higher_groups")) |>
  mutate(var_group=case_when(str_detect(variable, "pvr") ~ "phylogenetic",
                             str_detect(variable, "hfp") ~ "human footprint",
                             str_detect(variable, "biome") ~ "biome",
                             TRUE ~ variable)) |>
  group_by(var_group, fold) |>
  summarise(mean_decrease_accuracy=sum(mean_decrease_accuracy), .groups="drop") |>
  mutate(var_group=recode(var_group, "L3_count"="no. regions",
                          "humphreys_lifeform"="lifeform", "year"="year of description")) |>
  mutate(var_group=str_to_sentence(var_group)) |>
  group_by(var_group) |>
  median_qi(mean_decrease_accuracy, .width=0.95) |>
  mutate(var_group=reorder(var_group, mean_decrease_accuracy)) |>
    ggplot(mapping=aes(x=mean_decrease_accuracy, y=var_group, xmin=.lower, xmax=.upper, colour=mean_decrease_accuracy)) +
    geom_pointinterval(interval_colour="grey80", show.legend=FALSE,
                       fatten_point=2) +
    geom_vline(xintercept=0, linetype=2, size=1, color="red") +
    labs(x="Mean decrease in accuracy", y="") +
  viridis::scale_color_viridis(option="A", end=0.9) +
  theme(panel.background=element_blank(),
        panel.grid.major.x=element_line(linetype=2, colour="grey80"),
        panel.grid.major.y=element_blank(),
        text=element_text(family="DejaVu Sans"))

ggsave("output/bart-feature-importance-grouped.svg", grouped_importance_plot, height=3, width=6)

d <- tibble(L3_count=1:50)
t <- test_preds |>
  filter(cv == "Random CV") |>
  mutate(correct=as.numeric(obs == .pred_class)) |>
  filter(obs == "threatened")

fit <- glm(correct ~ log(L3_count), family=binomial(link="logit"), data=t)

line <- augment(fit, newdata=d, type.predict="response")

n_given_p <- function(p) {
  logit_p <- log(p / (1 - p))
  log_a <- (logit_p - fit$coefficients[[1]]) / fit$coefficients[[2]]
  exp(log_a)
}

glue::glue("Probability of correctly classifying a species with 1 region: {fmt(line[1,]$.fitted)}")
glue::glue("50% of threatened species with {floor(n_given_p(0.5))} regions are correctly classified.")
glue::glue("10% of threatened species with {floor(n_given_p(0.9))} regions are incorrectly classified.")

threatened_size <-
  distributions |>
    tibble() |>
    inner_join(pred_data |> filter(obs == "threatened") |> select(plant_name_id, L3_count), by="plant_name_id") |>
  group_by(area_code_l3) |>
  summarise(`Median range size (no. regions)`=median(L3_count), n_threat=n())

threat_data <- 
  country_perf |>
  left_join(threatened_size, by="area_code_l3") |>
  filter(cv == "Random CV", .metric == "sensitivity")

binom_model <- glm(.estimate ~ log(`Median range size (no. regions)`), family=binomial(link="logit"),
                   weights=n_threat, data=threat_data)

binom_line <- augment(binom_model, newdata=rename(d, `Median range size (no. regions)`=L3_count), type.predict="response")

sensitivity_vs_range <-
threat_data |>
  ggplot(mapping=aes(x=`Median range size (no. regions)`, y=.estimate)) +
  geom_point(alpha=0.5, colour="grey50") +
  geom_line(data=binom_line, mapping=aes(y=.fitted), size=1.5,
            colour="#ffae49") +
  labs(x="Median range size of a threatened species",
       y="Sensitivity") +
  theme(panel.background=element_blank(),
        panel.grid.major.x=element_line(linetype=2, colour="grey80"),
        panel.grid.major.y=element_line(linetype=2, colour="grey80"),
        text=element_text(family="DejaVu Sans")
        )

ggsave("output/bart-senitivity-range.svg", sensitivity_vs_range, height=3, width=6)

t <- test_preds |>
  filter(cv == "Random CV") |>
  mutate(correct=as.numeric(obs == .pred_class)) |>
  filter(obs == "threatened")
d <- tibble(area=seq(min(t$area), max(t$area), length.out=101))

fit <- glm(correct ~ log10(area), family=binomial(link="logit"), data=t)

line <- augment(fit, newdata=d, type.predict="response")

area_given_p <- function(p) {
  logit_p <- log(p / (1 - p))
  log_a <- (logit_p - fit$coefficients[[1]]) / fit$coefficients[[2]]
  10^(log_a)
}

glue::glue("Probability of correctly classifying a species with 1 smallest region: {fmt(predict(fit, newdata=tibble(area=min(t$area)), type='response')[[1]])}")
glue::glue("Probability of correctly classifying a species with 1 average region: {fmt(predict(fit, newdata=tibble(area=median(t$area)), type='response')[[1]])}")
glue::glue("Probability of correctly classifying a species with 1 large region: {fmt(predict(fit, newdata=tibble(area=max(t$area)), type='response')[[1]])}")
glue::glue("50% of threatened species with {floor(area_given_p(0.5))} km^2 range are correctly classified.")
glue::glue("10% of threatened species with {floor(area_given_p(0.9))} km^2 range are incorrectly classified.")

# criteria ----
full_rl <- read_csv("01_raw_data/redlist-all-vars-Jul2022_wcvpNewPhyt.csv", show_col_types=FALSE, progress=FALSE)

criteria <-
  full_rl |>
  select(accepted_plant_name_id, criteria)

threatened_taxa <-
  test_preds |>
  filter(cv == "Random CV",
         obs == "threatened") |>
  left_join(criteria, by=c("plant_name_id"="accepted_plant_name_id"))

split_criteria <- function(x) {
  cleaned <- str_remove_all(x, "[a-z\\,\\(\\)]+")
  str_split(cleaned, "[\\;\\, ]+(?=A|B|C|D)")
}

expand_criteria <- function(x) {
  first <- str_extract(x, "^[A-E]")
  split <- str_remove(x, "^[A-E]")
  split <- str_split(split, "\\+")
  split <- purrr::map(split, str_trim)
  purrr::map2(split, first, ~paste0(.y, .x))
}

threatened_taxa <-
  threatened_taxa |>
  filter(!is.na(criteria)) |>
  mutate(criteria_str=criteria,
    criteria=split_criteria(criteria),
         ) |>
  unnest(criteria) |>
  rowwise() |>
  mutate(criteria=expand_criteria(criteria)) |>
  unnest(criteria) |>
  filter(criteria != "", !is.na(criteria))

threatened_criteria <-
  threatened_taxa |>
  mutate(criteria=str_extract(criteria, "^[A-E]")) |>
  mutate(cited="cited") |>
  distinct(id, plant_name_id, obs, .pred_class, criteria, cited) |>
  pivot_wider(names_from=criteria, values_from=cited, 
              values_fill="not cited") |>
  pivot_longer(c(A, B, C, D), names_to="criterion", values_to="cited") 

cite_proportions <-
  threatened_criteria |>
    filter(cited == "cited") |>
    group_by(criterion) |>
    summarise(p=n() / length(unique(threatened_criteria$plant_name_id)), .groups="drop") |>
  mutate(label=glue::glue("{criterion} ({scales::label_percent(accuracy=0.1)(p)})"))

criteria_plot <- 
threatened_criteria |>
  disaggregate_performance(id, criterion, cited, metrics=accuracy) |>
  group_by(criterion, cited) |>
  median_qi(.estimate) |>
  left_join(cite_proportions, by=c("criterion")) |>
  mutate(cited=factor(cited, levels=c("not cited", "cited"))) |>
  ggplot(mapping=aes(x=.estimate, xmin=.lower, xmax=.upper, y=cited)) +
  geom_pointinterval(mapping=aes(colour=cited)) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_colour_manual(values=c(cited="#41B0AB", `not cited`="grey80")) +
  guides(alpha="none", colour="none") +
  facet_wrap(~label) +
  labs(x="Sensitivity", y="") +
  theme(panel.background=element_blank(),
        panel.border=element_rect(colour="grey80", fill=NA),
        panel.grid.major.x=element_line(linetype=2, colour="grey80"),
        panel.grid.major.y=element_blank(),
        text=element_text(family="DejaVu Sans"))

ggsave("output/bart-criteria-performance.svg", criteria_plot, height=3, width=6)

b_taxa <- 
  threatened_taxa |>
  filter(str_detect(criteria, "B")) |>
  distinct(id, plant_name_id, obs, .pred_class, criteria) |>
  mutate(cited=TRUE) |>
  pivot_wider(names_from=criteria, values_from=cited, values_fill=FALSE) |>
  mutate(b_criteria=case_when(B1 & B2 ~ "both",
                              B1 ~ "B1",
                              B2 ~ "B2"))
  

b_proportions <-
  b_taxa |>
    group_by(b_criteria) |>
    summarise(p=n() / nrow(b_taxa), .groups="drop") |>
  mutate(label=glue::glue("{b_criteria} ({scales::label_percent(accuracy=0.1)(p)})"))


b_criteria_plot <-
b_taxa |>
  disaggregate_performance(id, b_criteria) |>
  group_by(b_criteria) |>
  median_qi(.estimate) |>
  left_join(b_proportions, by="b_criteria") |>
  mutate(label=factor(label),
         label=fct_rev(label)) |>
  ggplot(mapping=aes(x=.estimate, xmin=.lower, xmax=.upper, 
                     y=label)) +
  geom_pointinterval() +
  geom_text(mapping=aes(label=label), hjust=1.25) +
  scale_x_continuous(limits=c(0, 1)) +
  labs(x="Sensitivity", y="") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.x=element_line(linetype=2, colour="grey80"),
        panel.grid.major.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(family="DejaVu Sans"))

ggsave("output/bart-b-criterion-performance.svg", b_criteria_plot, height=3, width=6)

# threats ----
threats <- read_csv("01_raw_data/threats.csv", show_col_types=FALSE)
threatened_threats <-
  full_rl |>
  select(taxonid, accepted_plant_name_id) |>
  inner_join(
    threats, 
    by=c("taxonid"="internalTaxonId")
  ) |>
  inner_join(
    threatened_taxa |>
      distinct(plant_name_id, obs, .pred_class),
    by=c("accepted_plant_name_id"="plant_name_id")
  ) |>
  mutate(correct=ifelse(obs == .pred_class, "correct", "incorrect"))

nclass <- 
  threatened_threats |>
  group_by(correct) |>
  summarise(n_species=n_distinct(taxonid))

threat_comparison <-
threatened_threats |>
  mutate(major_code=str_extract(code, "^\\d+")) |>
  distinct(taxonid, correct, major_code) |>
  count(correct, major_code) |>
  left_join(nclass, by="correct") |>
  mutate(p=n / n_species) |>
  mutate(major_code=as.numeric(major_code)) |>
  ggplot(mapping=aes(x=factor(major_code), y=p, fill=correct)) +
  geom_col(position="dodge") +
  scale_fill_manual(values=c("correct"="grey80", 
                             "incorrect"="#E77770"),
                    name="") +
  labs(x="Major threat code", y="Proportion of threatened species affected") +
  theme(legend.position="bottom",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.y=element_line(linetype=2, colour="grey80"),
        panel.grid.major.x=element_blank(),
        text=element_text(family="DejaVu Sans"))

ggsave("output/bart-threat-comparison.svg", threat_comparison, height=3, width=6)

# predictions ----

overall_preds <- 
  pred_data |>
  mutate(name="Overall") |>
  plot_threat(name, draws=predictions) +
  scico::scale_colour_scico(
    palette="lajolla",
    na.value="grey80",
    limits=c(0, 1)
  ) +
  remove_xaxis()

zone_preds <-
  pred_data |>
  left_join(species_zones,
            by="plant_name_id") |>
  plot_threat(zone, draws=predictions) +
  scico::scale_colour_scico(
    palette="lajolla",
    na.value="grey80",
    limits=c(0, 1)
  )

climate_preds <-
  pred_data |>
  replace_na(list(climate_description="Unknown")) |>
  plot_threat(climate_description, draws=predictions) +
  scico::scale_colour_scico(
    palette="lajolla",
    na.value="grey80",
    limits=c(0, 1)
  )

lifeform_preds <-
  pred_data |>
  replace_na(list(humphreys_lifeform="unknown")) |>
  mutate(humphreys_lifeform=str_to_sentence(humphreys_lifeform)) |>
  plot_threat(humphreys_lifeform, draws=predictions) +
  scico::scale_colour_scico(
    palette="lajolla",
    na.value="grey80",
    limits=c(0, 1)
  ) +
  remove_xaxis()

family_preds <-
  pred_data |>
  filter(family %in% show_families) |>
  plot_threat(family, draws=predictions) +
  scico::scale_colour_scico(
    palette="lajolla",
    na.value="grey80",
    limits=c(0, 1)
  )

pred_plot_climate <- ((overall_preds / lifeform_preds / climate_preds) + 
     plot_layout(heights=c(1, 5, 8)) | family_preds) +
  plot_layout(guides="collect",
              heights=c(1, 1)) &
    theme(legend.position="bottom",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.x=element_line(linetype=2, colour="grey80"),
        panel.grid.major.y=element_blank(),
        text=element_text(family="DejaVu Sans"))

ggsave("output/bart-predictions-disaggregated-climate.svg", pred_plot_climate, height=10, width=7)

pred_plot_zones <- ((overall_preds / lifeform_preds / zone_preds) + 
     plot_layout(heights=c(1, 5, 4)) | family_preds) +
  plot_layout(guides="collect",
              heights=c(1, 1)) &
    theme(legend.position="bottom",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.x=element_line(linetype=2, colour="grey80"),
        panel.grid.major.y=element_blank(),
        text=element_text(family="DejaVu Sans"))

ggsave("output/bart-predictions-disaggregated-zones.svg", pred_plot_zones, height=10, width=7)

country_preds <-
    distributions |>
    left_join(
      predictions |> select(plant_name_id, .draw, threatened),
      by="plant_name_id"
    ) |>
    group_by(area_code_l3, .draw)

country_preds <-
  country_preds |>
  summarise(threatened=mean(threatened, na.rm=TRUE)) |>
  median_qi(threatened, na.rm=TRUE)

threat_map <-
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
  )

uncertain_map <-
  country_preds |>
  mutate(uncertainty=.upper - .lower) |>
  plot_map(uncertainty, .proj="moll", .points=TRUE) +
  scico::scale_fill_scico(
    palette="grayC",
    na.value="grey80",
    labels=scales::label_percent(),
    name="Uncertainty (width of 95 % CI)"
  ) +
  scico::scale_colour_scico(
    palette="grayC",
    na.value="grey80",
    labels=scales::label_percent(),
    name="Uncertainty (proportion unresolved)"
  )

pred_map <- threat_map + uncertain_map & plot_annotation(tag_level="A") & theme(text=element_text(family="DejaVu Sans"))
ggsave("output/bart-prediction-map.svg", pred_map, height=5, width=10)

# spreadsheet ----

overall_perf_summary <- overall_performance |>
  mutate(blah=glue::glue("{fmt(.estimate)}\n[{fmt(.lower)}, {fmt(.upper)}]")) |>
  select(cv, .metric, blah) |>
  pivot_wider(names_from=.metric, values_from=blah)

criteria_perf <-
threatened_criteria |>
  disaggregate_performance(id, criterion, cited, metrics=accuracy) |>
  group_by(criterion, cited) |>
  median_qi(.estimate) |>
  left_join(cite_proportions, by=c("criterion"))

b_criteria_perf <-
b_taxa |>
  disaggregate_performance(id, b_criteria) |>
  group_by(b_criteria) |>
  median_qi(.estimate) |>
  left_join(b_proportions, by="b_criteria")

full_importance_summary <- importance |>
    group_by(variable) |>
    median_qi(.width=0.95) |>
    mutate(variable=reorder(variable, mean_decrease_accuracy)) 

importance_summary <- importance |>
  mutate(fold=c(rep(1, 88), rep(2, 88), 
                rep(3, 88), rep(4, 88), rep(5, 88))) |>
  filter(! variable %in% c("genus", "family", "order", "higher_groups")) |>
  mutate(var_group=case_when(str_detect(variable, "pvr") ~ "phylogenetic",
                             str_detect(variable, "hfp") ~ "footprint",
                             str_detect(variable, "biome") ~ "biome",
                             TRUE ~ variable)) |>
  group_by(var_group, fold) |>
  summarise(mean_decrease_accuracy=sum(mean_decrease_accuracy), .groups="drop") |>
  mutate(var_group=recode(var_group, "L3_count"="no. regions",
                          "humphreys_lifeform"="lifeform", "year"="year of description")) |>
  mutate(var_group=str_to_sentence(var_group)) |>
  group_by(var_group) |>
  median_qi(mean_decrease_accuracy, .width=0.95)

overall_summary <-
  predictions |> 
  left_join(
    pred_data |> mutate(name="Overall") |> select(plant_name_id, set, name), 
    by=c("plant_name_id", "set")
  ) |>
  group_by(name, .draw) |> 
  summarise(
    threatened=mean(threatened),
    .groups="drop"
  ) |>
  group_by(name) |>
  median_qi(threatened)

lifeform_summary <-
  predictions |> 
  left_join(
    pred_data |> 
    replace_na(list(humphreys_lifeform="unknown")) |>
    mutate(humphreys_lifeform=str_to_sentence(humphreys_lifeform)) |> 
    select(plant_name_id, set, humphreys_lifeform), 
    by=c("plant_name_id", "set")
  ) |>
  group_by(humphreys_lifeform, .draw) |> 
  summarise(
    threatened=mean(threatened),
    .groups="drop"
  ) |>
  group_by(humphreys_lifeform) |>
  median_qi(threatened)

family_summary <-
  predictions |> 
  left_join(
    pred_data |> select(plant_name_id, set, family), 
    by=c("plant_name_id", "set")
  ) |>
  group_by(family, .draw) |> 
  summarise(
    threatened=mean(threatened),
    .groups="drop"
  ) |>
  group_by(family) |>
  median_qi(threatened)

zone_data <-
  pred_data |>
  left_join(species_zones, by="plant_name_id")

zone_summary <-
  predictions |>
  left_join(
    zone_data |>  select(plant_name_id, set, zone), 
    by=c("plant_name_id", "set")
  ) |>
  group_by(zone, .draw) |> 
  summarise(
    threatened=mean(threatened),
    .groups="drop"
  ) |>
  group_by(zone) |>
  median_qi(threatened)

climate_summary <-
  predictions |> 
  left_join(
    pred_data |> 
    replace_na(list(climate_description="unknown")) |>
    select(plant_name_id, set, climate_description), 
    by=c("plant_name_id", "set")
  ) |>
  group_by(climate_description, .draw) |> 
  summarise(
    threatened=mean(threatened),
    .groups="drop"
  ) |>
  group_by(climate_description) |>
  median_qi(threatened)

write_xlsx(list(
  "perf-overall"=overall_perf_summary,
  "perf-lifeform"=filter(lifeform_perf, .metric=="TSS"),
  "perf-family"=filter(family_perf, .metric=="TSS"),
  "perf-zone"=filter(zone_perf, .metric=="TSS"),
  "perf-climate"=filter(climate_perf, .metric=="TSS"),
  "perf-crit"=criteria_perf,
  "perf-b-crit"=b_criteria_perf,
  "importance-full"=full_importance_summary,
  "importance-group"=importance_summary,
  "pred-overall"=overall_summary,
  "pred-lifeform"=lifeform_summary,
  "pred-family"=family_summary,
  "pred-zone"=zone_summary,
  "pred-climate"=climate_summary
), "output/bart-summary.xlsx")

# summarise predicted probabilities ----
ncores <- parallelly::availableCores()
cluster <- new_cluster(ncores)
cluster_library(cluster, c("tidyverse", "ggdist"))

species_probs <-
  probabilities |>
  group_by(set, plant_name_id) |>
  partition(cluster) |>
  mean_hdi() |>
  collect()

write_csv(species_probs, "output/bart-species-probabilities.csv")
