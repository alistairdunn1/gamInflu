#' Step plot
#'
#' A plot of the Standardised indices as each explanatory variable is added to the model
#'
#' @param obj The influence object
#' @param panels Whether to use panels or not
#' @param ... Additional arguments passed to ggplot
#' @return A ggplot object or a list of ggplot objects
#' @export
step_plot <- function(obj, panels = TRUE, ...) {
  if (is.null(obj$indices)) {
    stop("Call calc() before plotting")
  }
  
  # Determine which columns to plot
  start_col <- 6
  cols <- start_col:ncol(obj$indices)
  
  if (panels) {
    # Create a separate plot for each step
    plots <- list()
    
    for (i in seq_along(cols)) {
      col <- cols[i]
      col_name <- names(obj$indices)[col]
      
      # Check if this is an interaction term
      is_interaction <- col_name %in% names(obj$interactions) || 
                        (substr(col_name, 1, 1) == "+" && substr(col_name, 2, nchar(col_name)) %in% names(obj$interactions))
      
      # Get the term name without the "+" prefix
      term_name <- if (substr(col_name, 1, 1) == "+") substr(col_name, 2, nchar(col_name)) else col_name
      
      # Create a prettier label for the plot title
      if (is_interaction && term_name %in% names(obj$interactions)) {
        interaction_vars <- obj$interactions[[term_name]]$vars
        if (obj$interactions[[term_name]]$type %in% c("tensor.smooth", "tp", "ti")) {
          # For tensor smooths, use a nice format
          plot_title <- paste0("ti(", paste(interaction_vars, collapse=","), ")")
        } else {
          # For parametric interactions
          plot_title <- paste(interaction_vars, collapse=":")
        }
      } else {
        plot_title <- col_name
      }
      
      # Prepare data for this step
      plot_data <- obj$indices %>%
        dplyr::select(level, all_of(names(obj$indices)[cols[1]:col])) %>%
        mutate(level_num = as.integer(level)) %>%
        pivot_longer(
          cols = -c(level, level_num),
          names_to = "step",
          values_to = "value"
        )
      
      # Define color and linetype based on step relative to current step
      plot_data <- plot_data %>%
        mutate(
          step_num = match(step, names(obj$indices)[cols[1]:col]),
          step_type = case_when(
            step == col_name ~ "current",
            step == names(obj$indices)[col-1] & col > cols[1] ~ "previous",
            TRUE ~ "other"
          )
        )
      
      # Create the plot
      p <- ggplot(plot_data, aes(x = level_num, y = value, group = step)) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
        geom_line(aes(linetype = step_type, color = step_type)) +
        geom_point(data = filter(plot_data, step_type == "current"), size = 3, color = obj$colour) +
        scale_linetype_manual(values = c("current" = "solid", "previous" = "dashed", "other" = "dotted")) +
        scale_color_manual(values = c("current" = obj$colour, "previous" = obj$colour, "other" = "grey70")) +
        scale_x_continuous(breaks = plot_data$level_num, labels = plot_data$level) +
        labs(
          x = if (i == length(cols)) obj$labels[[obj$focus]] else "",
          y = "Index",
          title = plot_title
        ) +
        theme_minimal() +
        theme(
          legend.position = "none",
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
    # Create a single plot with all steps
    plot_data <- obj$indices %>%
      dplyr::select(level, all_of(names(obj$indices)[cols])) %>%
      mutate(level_num = as.integer(level)) %>%
      pivot_longer(
        cols = -c(level, level_num),
        names_to = "step",
        values_to = "value"
      )
    
    # Create better labels for steps
    plot_data <- plot_data %>%
      mutate(
        step_label = sapply(step, function(s) {
          term_name <- if (substr(s, 1, 1) == "+") substr(s, 2, nchar(s)) else s
          if (term_name %in% names(obj$interactions)) {
            interaction_vars <- obj$interactions[[term_name]]$vars
            if (obj$interactions[[term_name]]$type %in% c("tensor.smooth", "tp", "ti")) {
              paste0("ti(", paste(interaction_vars, collapse=","), ")")
            } else {
              paste(interaction_vars, collapse=":")
            }
          } else {
            s
          }
        })
      )
    
    # Create the plot
    ggplot(plot_data, aes(x = level_num, y = value, group = step_label, color = step_label)) +
      geom_line() +
      geom_point(size = 3) +
      scale_x_continuous(breaks = plot_data$level_num, labels = plot_data$level) +
      scale_color_brewer(palette = "Set1") +
      labs(
        x = obj$labels[[obj$focus]],
        y = "Index",
        color = "Step"
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      )
  }
}
