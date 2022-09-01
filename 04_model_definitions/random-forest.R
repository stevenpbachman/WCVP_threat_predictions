# MODEL SPECIFICATION ----
specify_model <- function() {
  rand_forest(
    trees=1000,
    mtry=tune(),
    min_n=tune()
  ) |>
    set_mode("classification") |>
    set_engine("randomForest")
}

# PRE-PROCESSING SPECIFICATION ----
specify_recipe <- function(data, ...) {
  form <- formula(
    obs ~ L3_count + humphreys_lifeform + climate_description
  )
  rec <- recipe(form, data=data) |>
    step_log(L3_count) |>
    step_unknown(humphreys_lifeform)
  
  rec
}

# HYPERPARAMETERS ----
hparam_grid <- 
  grid_regular(
    min_n(),
    mtry(range=c(1, 3)),
    levels=3
  )