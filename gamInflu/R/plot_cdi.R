#' CDI Plot for GAM Terms
#'
#' Create a Coefficient-Distribution-Influence plot for a smooth term's
#' influence on the parametric focus term. The plot consists of three panels:
#' \itemize{
#'   \item Top panel: Coefficients/effects of the plotted term
#'   \item Bottom left: Distribution of observations across term levels vs focus levels
#'   \item Bottom right: Standardized influence scores (mean = 1)
#' }
#'
#' @param x A gam_influence object with calculations performed
#' @param term The smooth term to plot (from expanded terms)
#' @param variable Optional variable name if it can't be inferred from term
#' @param type Optional type of plot (e.g., "coefficients", "distribution", "influence", or "all")
#' @param ... Additional plotting arguments
#' @return A combined ggplot object
#' @export
#' 
plot_cdi <- function(x, term, variable = NULL, type = "all", ...) {
  
  if (!x$calculated) {
    stop("Must run calculate() first")
  }
   
  if (!term %in% names(x$expanded_terms)) {
    stop("Term '", term, "' not found in expanded terms. Available terms: ",
         paste(names(x$expanded_terms), collapse = ", "))
  }
  
  # Get term information
  term_info <- x$expanded_terms[[term]]
  
  # Determine variable for distribution plot
  if (is.null(variable)) {
    variable <- infer_variable_for_cdi(x, term_info)
  }
  
  # Create individual plots
  valid_type <- c("coefficients", "distribution", "influence", "all")
  type <- valid_type[pmatch(type, valid_type, nomatch = 0)]
  if (!type %in% c("coefficients", "distribution", "influence", "all")) {
    stop("Invalid type specified. Choose from 'coefficients', 'distribution', 'influence', or 'all'.")
  } 

  if(type == "coefficients") {
    coef_plot <- plot_cdi_coefficients(x, term, ...)
    return(coef_plot) 
  } else if(type == "distribution") {
    dist_plot <- plot_cdi_distribution(x, variable, ...)
    return(dist_plot) 
  } else if(type == "influence") {
    infl_plot <- plot_cdi_influence(x, term, ...)
    return(infl_plot)
  } else if(type == "all") {
    coef_plot <- plot_cdi_coefficients(x, term, ...)
    dist_plot <- plot_cdi_distribution(x, variable, ...)
    infl_plot <- plot_cdi_influence(x, term, ...)

    combined_plot <- grid.arrange(
      coef_plot, arrangeGrob(dist_plot, infl_plot, ncol = 2),
      ncol = 1, heights = c(1, 1.2))
    return(combined_plot) 
  }
  stop("No valid plot type specified. Choose from 'coefficients', 'distribution', 'influence', or 'all'.")
}
