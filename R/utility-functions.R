#' Check for uncommited changes and ask user to continue.
#' 
#' @param .stop Whether to raise an error for uncommitted changes or just a warning.
#'   If in an interactive session, `.stop=TRUE` will give the user the option to continue
#'   or not.
#'
warn_uncommitted <- function(.stop=TRUE) {
  status <- git2r::status()
  uncommitted <- (length(status$staged) + length(status$unstaged)) > 0
  
  msg <- c(
    "{.strong Uncommited changes in your repository.}\n",
    "Commiting your changes before training your model ensures reproducibility."
  )
  
  if (.stop & !rlang::is_interactive() & uncommitted) {
    cli::cli_abort(msg)
  }
  
  if (!.stop & uncommitted) {
    cli::cli_warn(msg)
  }
  
  if (.stop & uncommitted)  {
    cli::cli_alert_warning(c(msg, "\nContinue without committing your changes?"))
    qs <- c("Yes", "Not now", "No way")
    qs <- sample(qs)
    
    out <- utils::menu(qs)
    .stop <- qs[[out]] != "Yes"
  }
  
  if (.stop & uncommitted) {
    cli::cli_abort("Stopping so user can commit changes.")
  }
  
  invisible()
}

#' Calculate test set performance across subgroups.
#' 
#' @param d A data frame of test set predictions.
#' @param ... Columns containing the names of the subgroups to calculate performance
#'  across.
#' @param metrics A metric or set of metrics from `yardstick` for evaluating performance.
#' 
#' @return A data frame
#'
disaggregate_performance <- function(d, ..., metrics=NULL) {
  if (is.null(metrics)) {
    metrics <- accuracy
  }
  
  join_keys <- sapply(enquos(...), rlang::quo_name)
  
  n <- count(d, ...)
  d |>
    group_by(...) |>
    metrics(truth=obs, estimate=.pred_class, event_level="second") |>
    left_join(n, by=join_keys)
}

#' Extract observation ids for training and test sets of a cross-validation fold.
#'
#' @param split an `rsample::rsplit` object for one cross-validation fold.
#'
#' @return a list with the plant_name_id for species in each set.
#'
get_fold_ids_ <- function(split) {
  train <- analysis(split)
  test <- assessment(split)

  list(
    train=train$plant_name_id,
    test=test$plant_name_id
  )
}

#' Extract observation ids for the training and test sets of all folds in a CV scheme.
#'
#' @param splits an `rsample::rset` object for a cross-validation scheme.
#'
#' @return a list with a list for each fold with the plant_name_id for species in the train and test sets.
#'
get_cv_ids <- function(splits) {
  ids <- lapply(splits$splits, get_fold_ids_)
  for (i in seq_along(splits$id)) {
    ids[[i]][["name"]] <- splits$id[[i]]
  }

  ids
}

#' Insert a key-value pair in every list in a list of named lists.
#'
#' @param lists a list of named lists.
#' @param items a named list of items to insert into every named list in `lists`.
#'
#' @ return a list of named lists.
#'
insert_items <- function(lists, items) {
  for (i in seq_along(lists)) {
    lists[[i]] <- c(lists[[i]], items)
  }

  lists
}