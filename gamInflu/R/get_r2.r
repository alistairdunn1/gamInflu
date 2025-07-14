#' Get Model r2
#'
#' @param x gam_influence object
#' @return The model_progression object from the calculated gam_influence object
#' @export
#'
get_r2 <- function(x) {
  if (is.null(x$model_progression)) {
    stop("Model progression has not been calculated. Run calculate() first.")
  }
  return(x$model_progression)
}
