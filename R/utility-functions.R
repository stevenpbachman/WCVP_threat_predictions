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
