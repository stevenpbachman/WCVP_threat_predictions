#' Functions to help train and evaluate the model

set_model_engine("bart", "classification", eng="mbart")
set_dependency("bart", eng="mbart", pkg="BART")

set_model_arg(
  model="bart",
  eng="mbart",
  parsnip="trees",
  original="ntree",
  func=list(pkg="dials", fun="trees", range=c(25, 250)),
  has_submodel=FALSE
)
set_model_arg(
  model="bart",
  eng="mbart",
  parsnip="prior_terminal_node_coef",
  original="base",
  func=list(pkg="dials", fun="prior_terminal_node_coef"),
  has_submodel=FALSE
)
set_model_arg(
  model="bart",
  eng="mbart",
  parsnip="prior_terminal_node_expo",
  original="power",
  func=list(pkg="dials", fun="prior_terminal_node_expo"),
  has_submodel=FALSE
)
set_model_arg(
  model="bart",
  eng="mbart",
  parsnip="prior_outcome_range",
  original="k",
  func=list(pkg="dials", fun="prior_outcome_range"),
  has_submodel=FALSE
)

set_fit(
  model="bart",
  eng="mbart",
  mode="classification",
  value=list(
    interface="data.frame",
    protect=c("x", "y"),
    func=c(fun="mbart_train"),
    defaults=list()
  )
)

set_encoding(
  model="bart",
  eng="mbart",
  mode="classification",
  options=list(
    predictor_indicators="one_hot",
    compute_intercept=FALSE,
    remove_intercept=TRUE,
    allow_sparse_x=FALSE
  )
)

set_pred(
  model="bart",
  eng="mbart",
  mode="classification",
  type="class",
  value=list(
    pre=NULL,
    post=NULL,
    func=c(fun="mbart_predict_calc"),
    args=list(
      obj=quote(object),
      new_data=quote(new_data),
      type="class"
    )
  )
)

set_pred(
  model="bart",
  eng="mbart",
  mode="classification",
  type="prob",
  value=list(
    pre=NULL,
    post=NULL,
    func=c(fun="mbart_predict_calc"),
    args=list(
      obj=quote(object),
      new_data=quote(new_data),
      type="prob"
    )
  )
)

set_pred(
  model="bart",
  eng="mbart",
  mode="classification",
  type="raw",
  value=list(
    pre=NULL,
    post=NULL,
    func=c(fun="mbart_predict_calc"),
    args=list(
      obj=quote(object),
      new_data=quote(new_data),
      type="ev"
    )
  )
)

set_pred(
  model="bart",
  eng="mbart",
  mode="classification",
  type="numeric",
  value=list(
    pre=NULL,
    post=NULL,
    func=c(fun="mbart_predict_calc"),
    args=list(
      obj=quote(object),
      new_data=quote(new_data),
      type="ppd"
    )
  )
)


permutation_importance_ <- function(fit_obj, test_set, .threshold=0.5, .n=1, parallel=FALSE) {
  if ("_bart" %in% class(fit_obj)) {
    fcn <- function(object, newdata) {
      p <- colMeans(predict(object$fit, newdata))
      ifelse(p > .threshold, "threatened", "not threatened")
    }
  } else if (length(levels(test_set$obs)) == 2) {
    fcn <- function(object, newdata) {
      p <- predict(object, newdata, type="prob")
      ifelse(p > .threshold, "threatened", "not threatened")
    }
  } else {
    fcn <- function(object, newdata) {
      p <- predict(object, newdata)
      p$.pred_class
    }
  }

  vip::vi_permute(
    fit_obj,
    train=test_set,
    target="obs",
    metric="accuracy",
    pred_wrapper=fcn,
    nsim=.n,
    parallel=parallel
  )
}

#' Given predicted class probabilites, choose the classification threshold that maximises
#' the TSS.
#' 
#' @preds A vector of predicted class probabilities.
#' @labels A vector of observed class labels corresponding to the predictions.
#' 
#' @return A list holding the threshold values tried (`alpha`), the TSS, and the 
#'  chosen threshold (`best`).
#'  
choose_threshold <- function(preds, labels) {
  rocr_preds <- prediction(preds, labels, label.order=c("not threatened", "threatened"))
  rocr_perf <- performance(rocr_preds, "sens", "spec")
  
  tss <- rocr_perf@x.values[[1]] + rocr_perf@y.values[[1]] - 1
  alpha <- rocr_perf@alpha.values[[1]]
  list(
    alpha=alpha,
    tss=tss,
    best=alpha[which.max(tss)]
  )
}

extract_samples <- function(obj, newdata, type="ev", ids=NULL) {
  if ("_mbart" %in% class(obj)) {
    samples <- mbart_samples(obj, newdata=newdata, type=type, ids=ids)
  } else {
    samples <- dbart_samples(obj, newdata=newdata, type=type, ids=ids)
  }

  samples
}

dbart_samples <- function(obj, newdata, type="ev", ids=NULL) {
  fit_obj <- extract_fit_engine(obj)
  pred <- predict(fit_obj, newdata=newdata, type=type)

  samples <- pred |>
    t() |>
    as_tibble()

  if (! is.null(ids)) {
    samples$.id <- ids
  }

  samples
}

mbart_samples <- function(obj, newdata, type="ev", ids=NULL) {
  if (type == "ev") {
    type <- "raw"
  } else {
    type <- "pred_int"
  }
  samples <- predict(obj, newdata=newdata, type=type)

  if (! is.null(ids)) {
    id_vec <- ids[samples$.id]
    samples$.id <- id_vec
  }

  samples
}

mbart_train <- function(x, y, ntree=50, k=2, power=2, base=0.95, sparse=FALSE, printevery=100000L) {
  y <- as.numeric(y)

  BART::mbart(x, y, ntree=ntree, k=k, power=power, base=base, sparse=sparse, printevery=printevery)
}

mbart_predict_calc <- function(obj, new_data, type) {
  types <- c("prob", "class", "ev", "ppd")
  
  mbart_obj <- extract_fit_engine(obj)
  
  preds <- predict(mbart_obj, as.data.frame(new_data))
  mn <- preds$prob.test.mean
  
  idx <- ceiling(seq(1, length(mn)) / length(obj$lvl))
  probs <- split(mn, idx)
  
  if (type == "prob") {
    res <- do.call(rbind, probs)
    res <- as_tibble(res)
    colnames(res) <- paste0(".pred_", obj$lvl)
  } else if (type == "class") {
    lvl_idx <- sapply(probs, which.max)
    lvl <- obj$lvl[lvl_idx]
    lvl <- factor(lvl, levels=obj$lvl)
    res <- tibble(.pred_class=lvl)
  } else {
    probs_tbl <- 
    preds$prob.test |>
    t() |>
    as_tibble()

    res <- split(probs_tbl, idx)

    if (type == "ppd") {
      res <- lapply(
        res,
        function(x) lapply(x, function(xi) obj$lvl[which.max(rmultinom(1, 1, xi))])
      )

      res <- lapply(res, as_tibble)
    } else {
      res <- lapply(res, function(x) mutate(x, .cls=obj$lvl))
    }

    res <- purrr::map2(res, seq(1, max(idx)), ~mutate(.x, .id=.y))
    
    bind_rows(res)
  }
  
  res
}
