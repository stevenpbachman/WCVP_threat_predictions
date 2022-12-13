#' Calculate conformal scores from binary classification probabilities.
#' 
calculate_scores <- function(d, label_col=NULL) {
  probs <- select(d, label=label_col, starts_with(".pred_"), -.pred_class)
  scores <- mutate(probs, across(starts_with(".pred_"), ~1-.x))
  
  if (! is.null(label_col)) {
    scores <- 
      scores |>
      pivot_longer(-label, names_to=".pred", values_to=".score", names_prefix=".pred_") |>
      filter(label == .pred) |>
      select(.score)
  } else {
    scores <- rename_with(scores, ~str_replace(.x, ".pred", ".score"))  
  }
  
  bind_cols(d, scores)
}

#' Calculate the score threshold for the desired significance level.
#' 
calculate_threshold <- function(scores, alpha) {
  n <- length(scores)
  threshold <- ceiling((n + 1) * (1 - alpha)) / n
  quantile(scores, threshold)[[1]]
}

#' Generate a prediction set from a calibrated score threshold.
make_prediction_set <- function(d, threshold) {
  probs <- select(d, starts_with(".pred_"), -.pred_class)
  pred_set <- mutate(probs, across(everything(), ~.x >= (1 - threshold)))
  pred_set <- rename_with(pred_set, ~str_replace(.x, ".pred", ".set"))
  
  bind_cols(d, pred_set)
}

#' Estimate the coverage by calibrating and testing on repeated splits of the data.
estimate_coverage <- function(d, alpha, times=100, label_col="obs", cls_conditional=FALSE) {
  eval_splits <- mc_cv(d, prop=0.5, times=times)
  
  eval_splits |>
    rowwise() |>
    mutate(qhat=calculate_threshold(analysis(splits)$.score, alpha=alpha)) |>
    mutate(coverage=list(estimate_coverage_(splits, alpha=alpha, label_col=label_col, 
                                            cls_conditional=cls_conditional))) |>
    select(resample_id=id, qhat, coverage) |>
    unnest(coverage)
}

#' Estimate the coverage on a single split of the data.
#' 
estimate_coverage_ <- function(split, alpha, label_col="obs", cls_conditional=FALSE) {
  cal <- analysis(split)
  val <- assessment(split)
  
  qhat <- calculate_threshold(cal$.score, alpha=alpha)
  sets <- make_prediction_set(val, qhat)
    
  
  coverage <- 
    sets |>
    select(starts_with(".set"), label=label_col) |>
    pivot_longer(-label, names_to=".pred", values_to="covered", names_prefix=".set_") |>
    filter(label == .pred)
  
  if (cls_conditional) {
    coverage |>
      group_by(label) |>
      summarise(coverage=mean(covered))
  } else {
    coverage <-
      coverage |>
      summarise(coverage=mean(covered))
    
    coverage[[label_col]] <- "both"
    coverage
  }
  
}
