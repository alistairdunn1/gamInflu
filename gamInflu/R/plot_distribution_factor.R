plot_distribution_factor <- function(focus_data, smooth_data, focus_name, variable_name) {
  
  # Create contingency table
  smooth_data <- as.factor(smooth_data)
  
  # Create data frame for ggplot
  dist_data <- data.frame(
    smooth_var = smooth_data,
    focus_var = focus_data
  )
  
  # Calculate proportions
  dist_summary <- dist_data %>%
    count(smooth_var, focus_var) %>%
    group_by(focus_var) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  p <- ggplot(dist_summary, aes(x = smooth_var, y = focus_var, size = prop)) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(name = "", guide = "none") +
    labs(x = variable_name, y = focus_name)
  
  return(p)
}
