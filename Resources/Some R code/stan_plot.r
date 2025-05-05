#' Standardization plot
#'
#' This plot compares the Standardised and unStandardised indices from a model
#'
#' @param obj The influence object
#' @param ... Additional arguments passed to ggplot
#' @return A ggplot object
#' @export
stan_plot <- function(obj, ...) {
  if (is.null(obj$indices)) {
    stop("Call calc() before plotting")
  }
  
  # Prepare data for plotting
  plot_data <- obj$indices %>%
    dplyr::select(level, unstan, stan, stanLower, stanUpper) %>%
    mutate(level_num = as.integer(level))
  
  # Create the plot
  ggplot(plot_data, aes(x = level_num)) +
    geom_errorbar(aes(ymin = stanLower, ymax = stanUpper), width = 0.2) +
    geom_line(aes(y = unstan, group = 1, color = "UnStandardised"), size = 1) +
    geom_point(aes(y = unstan, color = "UnStandardised"), size = 3) +
    geom_line(aes(y = stan, group = 1, color = "Standardised"), size = 1) +
    geom_point(aes(y = stan, color = "Standardised"), size = 3) +
    scale_x_continuous(breaks = plot_data$level_num, labels = plot_data$level) +
    scale_color_manual(values = c("UnStandardised" = "darkgrey", "Standardised" = obj$colour)) +
    labs(
      x = obj$labels[[obj$focus]],
      y = "Index",
      color = ""
    ) 
}
