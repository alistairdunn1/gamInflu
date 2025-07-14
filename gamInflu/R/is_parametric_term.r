#' Check if Term is Parametric
#'
#' @param model GAM model
#' @param term Term name to check
#' @return Logical indicating if term is parametric
#' @keywords internal
#' 
is_parametric_term <- function(model, term) {
  smooth_labels <- sapply(model$smooth, function(x) x$label)
  return(!term %in% smooth_labels)
}
