#' Influence Plot for GAM Terms
#'
#' Plot the influence of smooth terms on the parametric focus term
#'
#' @param x A gam_influence object with calculations performed
#' @param panels Whether to use separate panels for each term
#' @return A ggplot object
#' @export
plot_influence <- function(x, panels = TRUE) {
  if (!x$calculated) {
    stop("Must run calculate() first")
  }

  # Get influence data
  influence_data <- prepare_influence_data(x)

  if (panels) {
    return(plot_influence_panels(influence_data, x))
  } else {
    return(plot_influence_single(influence_data, x))
  }
}
