#' @title Subplot for Random Effects with By-Variable
#' @description Internal function to plot predicted effects for random effect terms with by-variables (e.g., s(vessel, by=f.NNets, bs="re")).
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param term_vars A character vector of variables in the term.
#' @param re_type Character; for random effects, one of "points", "qq", "hist", or "caterpillar".
#' @param cdi Logical indicating if the plot is for CDI (Cumulative Distribution Influence).
#' @return A ggplot object showing the random by-variable effects.
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar facet_wrap labs theme element_blank element_text xlab geom_hline geom_histogram geom_density geom_vline after_stat geom_abline
#' @importFrom patchwork plot_spacer
#' @importFrom stats coef vcov qnorm ppoints
#' @noRd
subplot_random_by_effect <- function(obj, t, term_vars, re_type = "points", cdi = FALSE) {
  message("Plotting random by-variable effects for term: ", t)

  # Extract by-variable from the term
  by_var <- sub(".*by\\s*=\\s*([^,\\)]+).*", "\\1", t)
  by_var <- trimws(gsub("[\"']", "", by_var))

  # Get the main variable (first in term_vars, typically the random effect variable)
  main_var <- term_vars[1]

  # Get model object to extract random effects
  model <- obj$model

  # Extract random effect coefficients for this term
  # For s(vessel, by=f.NNets, bs="re"), we need to find the relevant coefficients

  # Get all coefficient names that match this term pattern
  coef_names <- names(coef(model))

  # Find coefficients that match the random effect pattern for this term
  # Pattern for s(vessel, by=f.NNets, bs="re"): s(vessel):f.NNetsNet1.1, s(vessel):f.NNetsNet2.1, etc.
  pattern_base <- paste0("s\\(", main_var, "\\):", gsub("\\.", "\\\\.", by_var))
  matching_coefs <- grep(pattern_base, coef_names, value = TRUE)

  if (length(matching_coefs) == 0) {
    message("No matching random effect coefficients found for term: ", t)
    message("Looking for pattern: ", pattern_base)
    message("Available coefficient names: ", paste(head(coef_names, 10), collapse = ", "))
    return(patchwork::plot_spacer())
  }

  # Extract coefficients and their standard errors
  coef_values <- coef(model)[matching_coefs]
  coef_se <- sqrt(diag(vcov(model)))[matching_coefs]

  # Parse the coefficient names to extract levels and groups
  # Extract the by-variable level from coefficient names
  # by_levels <- sub(paste0(".*", gsub("\\.", "\\\\.", by_var)), "", matching_coefs)

  # Extract the main variable level (vessel level)
  # main_levels <- sub(paste0("s\\(", main_var, "\\)\\..*"), "", matching_coefs)
  # main_levels <- sub(paste0("s\\(", main_var, "\\)"), "", matching_coefs)
  # main_levels <- sub(paste0("\\.", gsub("\\.", "\\\\.", by_var), ".*"), "", main_levels)

  # Create a cleaner approach - extract levels from coefficient names
  if (main_var %in% names(obj$data) && by_var %in% names(obj$data)) {
    by_var_levels <- levels(as.factor(obj$data[[by_var]]))

    # Create a data frame for plotting
    plot_data <- data.frame()

    for (by_level in by_var_levels) {
      # Pattern: s(vessel):f.NNetsNet1.1, s(vessel):f.NNetsNet1.2, etc.
      by_pattern <- paste0("s\\(", main_var, "\\):", gsub("\\.", "\\\\.", by_var), by_level)
      by_coefs <- grep(by_pattern, matching_coefs, value = TRUE)

      if (length(by_coefs) > 0) {
        by_coef_values <- coef_values[by_coefs]
        by_coef_se <- coef_se[by_coefs]

        # Extract main variable indices from coefficient names
        # From "s(vessel):f.NNetsNet1.1" extract "1", from "s(vessel):f.NNetsNet1.2" extract "2", etc.
        main_indices <- sub(paste0(".*", gsub("\\.", "\\\\.", by_var), by_level, "\\."), "", by_coefs)

        # Create meaningful labels for main variable levels
        # Use the actual factor levels if available
        if (is.factor(obj$data[[main_var]])) {
          main_var_levels <- levels(obj$data[[main_var]])
          if (length(main_var_levels) >= max(as.numeric(main_indices), na.rm = TRUE)) {
            main_labels <- main_var_levels[as.numeric(main_indices)]
          } else {
            main_labels <- paste0(main_var, main_indices)
          }
        } else {
          main_labels <- paste0(main_var, main_indices)
        }

        temp_data <- data.frame(
          main_var = main_labels, # Use main_labels instead of main_levels_this
          by_var = by_level,
          coefficient = as.numeric(by_coef_values),
          se = as.numeric(by_coef_se),
          lower = as.numeric(by_coef_values) - 1.96 * as.numeric(by_coef_se),
          upper = as.numeric(by_coef_values) + 1.96 * as.numeric(by_coef_se),
          stringsAsFactors = FALSE
        )

        plot_data <- rbind(plot_data, temp_data)
      }
    }
  } else {
    # Fallback: create basic structure from coefficient names
    plot_data <- data.frame(
      main_var = seq_along(coef_values),
      by_var = "combined",
      coefficient = as.numeric(coef_values),
      se = as.numeric(coef_se),
      lower = as.numeric(coef_values) - 1.96 * as.numeric(coef_se),
      upper = as.numeric(coef_values) + 1.96 * as.numeric(coef_se),
      stringsAsFactors = FALSE
    )
  }

  if (nrow(plot_data) == 0) {
    message("No data available for plotting random by-variable effects")
    return(patchwork::plot_spacer())
  }

  # Create the plot based on re_type
  re_type <- pmatch(re_type, c("points", "qq", "hist", "caterpillar"), nomatch = 1)

  if (re_type == 1 || re_type == 4) { # points or caterpillar
    p <- ggplot(plot_data, aes(x = .data$main_var, y = .data$coefficient)) +
      geom_point(size = 2, colour = "royalblue") +
      geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper),
        width = 0.2, colour = "royalblue"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      facet_wrap(~by_var, scales = "free_x") +
      labs(
        y = "Random Effect Coefficient",
        title = paste("Random Effects:", main_var, "by", by_var)
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    if (re_type == 4) { # caterpillar - sort by coefficient value
      plot_data$main_var <- factor(plot_data$main_var,
        levels = plot_data$main_var[order(plot_data$coefficient)]
      )
      p <- p + aes(x = .data$main_var)
    }
  } else if (re_type == 2) { # qq plot - use simpler approach
    # Create Q-Q plot manually since stat_qq import might be problematic
    qq_data <- plot_data
    for (level in unique(plot_data$by_var)) {
      level_data <- plot_data[plot_data$by_var == level, ]
      n <- nrow(level_data)
      if (n > 1) {
        qq_data[qq_data$by_var == level, "theoretical"] <- qnorm(ppoints(n))
        qq_data[qq_data$by_var == level, "sample"] <- sort(level_data$coefficient)
      }
    }

    p <- ggplot(qq_data, aes(x = .data$theoretical, y = .data$sample)) +
      geom_point(colour = "royalblue") +
      geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
      facet_wrap(~by_var) +
      labs(
        x = "Theoretical Quantiles", y = "Sample Quantiles",
        title = paste("Q-Q Plot: Random Effects", main_var, "by", by_var)
      )
  } else if (re_type == 3) { # histogram
    p <- ggplot(plot_data, aes(x = .data$coefficient)) +
      geom_histogram(aes(y = after_stat(density)), bins = 10, alpha = 0.7, fill = "royalblue") +
      geom_density(colour = "royalblue") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      facet_wrap(~by_var) +
      labs(
        x = "Random Effect Coefficient", y = "Density",
        title = paste("Distribution: Random Effects", main_var, "by", by_var)
      )
  }

  # Remove axis labels for CDI plots
  if (cdi) {
    p <- p +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank()
      )
  } else {
    p <- p + xlab(main_var)
  }

  return(p)
}
