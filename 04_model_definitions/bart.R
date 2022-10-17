# MODEL SPECIFICATION ----
specify_model <- function() {
  bart(
    trees=tune(),
    prior_terminal_node_coef=tune(),
    prior_terminal_node_expo=tune(),
    prior_outcome_range=tune()
  ) |>
    set_mode("classification") |>
    set_engine("dbarts")
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
hparam_grid <- NULL