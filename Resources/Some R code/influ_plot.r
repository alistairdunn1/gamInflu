#' Influence plot
#'
#' A plot of the influence of each explanatory variable in the model
#'
#' @param obj The influence object
#' @param panels Whether to use panels or not
#' @param ... Additional arguments passed to ggplot
#' @return A ggplot object or a list of ggplot objects
#' @export
influ_plot <- function(obj, panels = TRUE, ...) {
  if (is.null(obj$influences)) {
    stop("Call calc() before plotting")
  }
  
  # Determine which columns to plot (skip the level column)
  cols <- 2:ncol(obj$influences)
  
  if (panels) {
    # Create a separate plot for each variable
    plots <- list()
    
    for (i in seq_along(cols)) {
      col <- cols[i]
      var_name <- names(obj$influences)[col]
      
      # Create a prettier label for the plot title
      if (var_name %in% names(obj$interactions)) {
        interaction_vars <- obj$interactions[[var_name]]$vars
        if (obj$interactions[[var_name]]$type %in% c("tensor.smooth", "tp", "ti")) {
          # For tensor smooths, use a nice format
          plot_title <- paste0("ti(", paste(interaction_vars, collapse=","), ")")
        } else {
          # For parametric interactions
          plot_title <- paste(interaction_vars, collapse=":")
        }
      } else {
        plot_title <- var_name
      }
      
      # Prepare data for this variable
      plot_data <- obj$influences %>%
        dplyr::select(level, all_of(var_name)) %>%
        mutate(
          level_num = as.integer(level),
          influence = exp(get(var_name))
        )
      
      # Create the plot
      p <- ggplot(plot_data, aes(x = level_num, y = influence)) +
        geom_hline(yintercept = 1, linetype = "dashed") +
        geom_line(color = obj$colour) +
        geom_point(size = 3, color = obj$colour) +
        scale_x_continuous(breaks = plot_data$level_num, labels = plot_data$level) +
        scale_y_log10() +
        labs(
          x = if (i == length(cols)) obj$labels[[obj$focus]] else "",
          y = "Influence",
          title = plot_title
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = if (i == length(cols)) 45 else 0, hjust = if (i == length(cols)) 1 else 0.5),
          axis.title.x = element_text(margin = margin(t = 10)),
          panel.grid.minor = element_blank()
        )
      
      # Add to list of plots
      plots[[i]] <- p
    }
    
    # Combine plots using patchwork
    combined_plot <- wrap_plots(plots, ncol = 1)
    return(combined_plot)
  } else {
    # Create a single plot with all variables
    plot_data <- obj$influences %>%
      dplyr::select(level, all_of(names(obj$influences)[cols])) %>%
      mutate(level_num = as.integer(level)) %>%
      pivot_longer(
        cols = -c(level, level_num),
        names_to = "variable",
        values_to = "log_influence"
      ) %>%
      mutate(influence = exp(log_influence))
    
    # Create better labels for variables
    plot_data <- plot_data %>%
      mutate(
        var_label = sapply(variable, function(v) {
          if (v %in% names(obj$interactions)) {
            interaction_vars <- obj$interactions[[v]]$vars
            if (obj$interactions[[v]]$type %in% c("tensor.smooth", "tp", "ti")) {
              paste0("ti(", paste(interaction_vars, collapse=","), ")")
            } else {
              paste(interaction_vars, collapse=":")
            }
          } else {
            v
          }
        })
      )
    
    # Create the plot
    ggplot(plot_data, aes(x = level_num, y = influence, group = var_label, color = var_label)) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_line() +
      geom_point(size = 3) +
      scale_x_continuous(breaks = plot_data$level_num, labels = plot_data$level) +
      scale_y_log10() +
      scale_color_brewer(palette = "Set1") +
      labs(
        x = obj$labels[[obj$focus]],
        y = "Influence",
        color = "Variable"
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      )
  }
}

