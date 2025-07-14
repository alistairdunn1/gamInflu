#' Plot Distribution for Numeric Smooth Variable (ggplot)
#'
#' @param focus_data Focus term data
#' @param smooth_data Smooth variable data
#' @param focus_name Focus term name
#' @param variable_name Variable name
#' @return ggplot object
#' @keywords internal
#' 
plot_distribution_numeric <- function(focus_data, smooth_data, focus_name, variable_name) {
  
  # Cut numeric variable into bins
  breaks <- pretty(smooth_data, n = 15)
  smooth_binned <- cut(smooth_data, breaks, include.lowest = TRUE, labels = FALSE)
  
  # Create data frame
  dist_data <- data.frame(
    smooth_var = smooth_binned,
    focus_var = focus_data
  )
  
  # Calculate proportions
  dist_summary <- dist_data %>%
    count(smooth_var, focus_var) %>%
    group_by(focus_var) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    filter(n > 0)
  
  # Extract bin midpoints for x-axis
  dist_summary$bin_mid <- breaks[as.numeric(dist_summary$smooth_var)]
  
  xx <<- dist_summary

  p <- ggplot(dist_summary, aes(x = bin_mid, y = focus_var, size = prop)) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(name = "") +
    labs(x = variable_name, y = focus_name)
  
  return(p)
}
