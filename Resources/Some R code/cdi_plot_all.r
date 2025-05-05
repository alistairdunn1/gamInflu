#' Plot a CDI plot for each of the terms in the model
#'
#' @param obj The influence object
#' @param ... Additional arguments passed to cdi_plot
#' @return A list of plots
#' @export
cdi_plot_all <- function(obj, ...) {
  plots <- list()
  
  for (term in obj$terms) {
    if (term != obj$focus) {
      plots[[term]] <- cdi_plot(obj, term, ...)
    }
  }
  
  return(plots)
}

