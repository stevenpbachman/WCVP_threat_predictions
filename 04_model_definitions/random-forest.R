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
  outcome <- "obs"
  predictors <- c(
    "L3_count",
    "humphreys_lifeform",
    "climate_description",
    "area",
    "year"
  )
  impute_phylo <- c(
    "higher_groups",
    "order",
    "family",
    "genus",
    "climate_description",
    "humphreys_lifeform"
  )

  rec <- 
    data |>
    select(all_of(c(outcome, predictors, impute_phylo)), starts_with("pvr"), starts_with("biome")) |>
    recipe() |>
    update_role(all_of(outcome), new_role="outcome") |>
    update_role(one_of(predictors), new_role="predictor") |>
    update_role(setdiff(impute_phylo, predictors), new_role="impute phylo") |>
    add_role(intersect(predictors, impute_phylo), new_role="impute phylo") |>
    step_unknown(
      humphreys_lifeform,
      climate_description
    ) |>
    step_impute_knn(
      starts_with("pvr"), 
      impute_with=imp_vars(has_role(match="impute phylo"))
    ) |>
    step_impute_bag(
      starts_with("biome"),
      impute_with=imp_vars(starts_with("pvr"))
    ) |>
    step_log(L3_count, area)
}

# HYPERPARAMETERS ----
make_grid <- function(data) {
  grid_regular(
    min_n(),
    finalize(mtry(), data),
    levels=3
  )
}
