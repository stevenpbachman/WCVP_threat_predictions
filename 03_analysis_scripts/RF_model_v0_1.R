
# following steps outlined in Baz's tidyassessments vignette
# https://barnabywalker.github.io/tidyassessments/articles/tidymodels-assessments.html

library(dplyr)
library(readr)
library(tidyr)
library(glue)
library(tidyverse)
library(scales)
library(tidymodels)
library(tidyassessments)

#### Input data ####
#split the species into ones that have been assessed (labelled) and those that haven’t (unlabelled). We’ll include Data Deficient species in unlabelled.
labelled <- filter(angio_pred_v1, category != "DD", ! is.na(category))
unlabelled <- filter(angio_pred_v1, category == "DD" | is.na(category))

glue("{nrow(labelled)} labelled examples for training and evaluation", 
     "{nrow(unlabelled)} unlabelled examples for prediction", .sep="\n")

# make our prediction target the binarised categories “threatened” (VU, EN, CR) and “not threatened” (LC and NT), to simplify things.
labelled$obs <- ifelse(labelled$category %in% c("LC", "NT"), "not threatened", "threatened")
labelled$obs <- factor(labelled$obs, levels=c("threatened", "not threatened"))

glue("{percent_format(accuracy=0.1)(mean(labelled$obs == 'threatened'))}",
     "of training/evaluation species are threatened", .sep=" ")

str(angio_pred_v1)

#### Data budget ####
folds <- vfold_cv(v=5, data=labelled)
folds

#### Evaluation metrics ####
metrics <- metric_set(accuracy, sensitivity, specificity, j_index)
metrics

#### Pre-processing ####
form <- formula(
  #obs ~ family + genus + humphreys_lifeform + climate_description + L3_count
  obs ~ L3_count + humphreys_lifeform + climate_description 
)

rec <- 
  recipe(form, data=labelled) #%>%
  #step_impute_knn(all_predictors()) %>%
  #step_corr(all_predictors(), threshold=0.9) %>%
  #step_zv(all_predictors()) %>%
  #step_normalize(all_predictors())

rec

#### Model specification ####
spec <-
  rand_forest(trees=1000) %>%
  set_engine("randomForest", importance=TRUE) %>%
  set_mode("classification")

spec

#### Workflow ####
wf <- 
  workflow() %>%
  add_recipe(rec) %>%
  add_model(spec)

wf

#### Model evaluation ####
results <- fit_resamples(
  wf, 
  folds, 
  metrics=metrics,
  control=control_resamples(save_pred=TRUE, save_workflow=TRUE)
)

collect_metrics(results)

#### Predict threat ####
m <- fit(wf, labelled)
m

preds <- predict(m, angio_pred_v1)
prop_threatened <- mean(preds$.pred_class == "threatened")
glue("{percent_format(accuracy=0.1)(prop_threatened)} of Angiosperm species predicted threatened")
