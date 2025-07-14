#' Calculate Fit Statistics
#'
#' @param x gam_influence object
#' @return List with fit statistics
#' @keywords internal
#' 
calculate_fit_statistics <- function(x) {
  model_summary <- summary(x$model)

  return(list(
    r_squared = model_summary$r.sq,
    adjusted_r_squared = model_summary$r.sq,
    deviance_explained = model_summary$dev.expl,
    aic = AIC(x$model),
    bic = BIC(x$model),
    gcv = model_summary$sp.criterion
  ))
}
