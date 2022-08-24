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
    step_log(L3_count)
  
  rec
}

# HYPERPARAMETERS ----
hparam_grid <- 
  grid_regular(
    trees(range=c(100, 300)),
    prior_terminal_node_coef(range=c(0.75, 0.95)),
    prior_terminal_node_expo(range=c(1, 3)),
    prior_outcome_range(range=c(1, 3)),
    levels=3
  )