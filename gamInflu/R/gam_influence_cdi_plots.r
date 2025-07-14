#' @title CDI plotting methods for GAMInfluence class
#' @description# CDI plot for 1D smooths
#' @param term The term to plot
#' @param variable The variable to plot (optional)
#' @param colour The colour to use for the plot
#' @param n_bins The number of bins to use for the plot
#' @param smooth_info Information about the smooth term
#' @return A list of ggplot objects for the CDI plots
#' @keywords internal
#'
cdi_plot_1d <- function(term, variable, colour, n_bins, smooth_info) {
  # Auto-detect variable if not provided
  if (is.null(variable)) {
    variable <- smooth_info$variables[1]
  }

  # Prepare coefficient data
  coeff_data <- private$prepare_coefficient_data_1d(term, variable, n_bins)

  # Prepare distribution data
  distrib_data <- private$prepare_distribution_data(term, variable, n_bins)

  # Prepare influence data
  if (term %in% names(self$influences)) {
    influence_data <- self$influences %>%
      dplyr::mutate(
        level_num = as.numeric(as.factor(level)),
        exp_influence = exp(!!rlang::sym(term))
      )
  } else {
    influence_data <- data.frame(
      level = self$influences$level,
      level_num = as.numeric(as.factor(self$influences$level)),
      exp_influence = rep(1, nrow(self$influences))
    )
  }

  # Create coefficient plot
  coeff_plot <- ggplot2::ggplot(coeff_data, ggplot2::aes(x = level_num, y = exp_coeff)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = exp_lower, ymax = exp_upper),
      fill = colour, alpha = 0.3
    ) +
    ggplot2::geom_line(colour = colour) +
    ggplot2::geom_point(colour = colour) +
    ggplot2::scale_x_continuous(
      breaks = coeff_data$level_num,
      labels = coeff_data$level
    ) +
    ggplot2::labs(
      title = paste("Smooth Effect:", smooth_info$description),
      x = private$get_label(variable),
      y = "Effect"
    )

  # Create distribution plot
  if (nrow(distrib_data) > 0) {
    distrib_plot <- ggplot2::ggplot(distrib_data, ggplot2::aes(x = term_level, y = focus_level)) +
      ggplot2::geom_point(ggplot2::aes(size = prop), colour = colour, alpha = 0.7) +
      ggplot2::scale_size_continuous(name = "Proportion", range = c(1, 8)) +
      ggplot2::labs(
        title = "Data Distribution",
        x = private$get_label(variable),
        y = private$get_label(self$focus)
      )
  } else {
    distrib_plot <- ggplot2::ggplot() +
      ggplot2::geom_text(ggplot2::aes(x = 0.5, y = 0.5, label = "No distribution data available"))
  }

  # Create influence plot
  influence_plot <- ggplot2::ggplot(influence_data, ggplot2::aes(x = exp_influence, y = level_num)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_line(colour = colour) +
    ggplot2::geom_point(colour = colour) +
    ggplot2::scale_y_continuous(
      breaks = influence_data$level_num,
      labels = influence_data$level
    ) +
    ggplot2::labs(
      title = "Influence on Focus Term",
      x = "Influence",
      y = private$get_label(self$focus)
    )

  return(list(
    coefficient = coeff_plot,
    distribution = distrib_plot,
    influence = influence_plot,
    combined = gridExtra::grid.arrange(coeff_plot, distrib_plot, influence_plot,
      layout_matrix = rbind(c(1, 1), c(2, 3))
    )
  ))
}

# CDI plot for 2D smooths
cdi_plot_2d <- function(term, variable, colour, n_bins, smooth_info) {
  variables <- smooth_info$variables

  if (length(variables) < 2) {
    stop("2D smooth requires at least 2 variables")
  }

  var1 <- variables[1]
  var2 <- variables[2]

  # Create prediction grid
  if (var1 %in% names(self$data) && var2 %in% names(self$data)) {
    var1_range <- range(self$data[[var1]], na.rm = TRUE)
    var2_range <- range(self$data[[var2]], na.rm = TRUE)

    pred_grid <- expand.grid(
      x1 = seq(var1_range[1], var1_range[2], length.out = 20),
      x2 = seq(var2_range[1], var2_range[2], length.out = 20)
    )
    names(pred_grid) <- c(var1, var2)

    # Add other variables at their means/modes
    for (var in names(self$data)) {
      if (!var %in% c(var1, var2, self$response)) {
        if (is.numeric(self$data[[var]])) {
          pred_grid[[var]] <- mean(self$data[[var]], na.rm = TRUE)
        } else {
          pred_grid[[var]] <- names(sort(table(self$data[[var]]), decreasing = TRUE))[1]
        }
      }
    }

    # Get predictions
    tryCatch(
      {
        preds <- predict(self$model, newdata = pred_grid, type = "terms", se.fit = TRUE)

        if (term %in% colnames(preds$fit)) {
          pred_grid$fit <- preds$fit[, term]
          pred_grid$se <- preds$se.fit[, term]
        } else {
          pred_grid$fit <- 0
          pred_grid$se <- 0
        }
      },
      error = function(e) {
        pred_grid$fit <- 0
        pred_grid$se <- 0
      }
    )

    # Create 2D effect plot
    effect_plot <- ggplot2::ggplot(pred_grid, ggplot2::aes_string(x = var1, y = var2, fill = "fit")) +
      ggplot2::geom_raster() +
      viridis::scale_fill_viridis_c(name = "Effect") +
      ggplot2::labs(
        title = paste("2D Smooth Effect:", smooth_info$description),
        x = private$get_label(var1),
        y = private$get_label(var2)
      )

    # Create contour plot
    contour_plot <- ggplot2::ggplot(pred_grid, ggplot2::aes_string(x = var1, y = var2, z = "fit")) +
      ggplot2::geom_contour_filled() +
      viridis::scale_fill_viridis_d(name = "Effect") +
      ggplot2::labs(
        title = "Effect Contours",
        x = private$get_label(var1),
        y = private$get_label(var2)
      )

    # Data distribution plot
    if (var1 %in% names(self$data) && var2 %in% names(self$data)) {
      data_plot <- ggplot2::ggplot(self$data, ggplot2::aes_string(x = var1, y = var2)) +
        ggplot2::geom_point(alpha = 0.5, colour = colour) +
        ggplot2::labs(
          title = "Data Distribution",
          x = private$get_label(var1),
          y = private$get_label(var2)
        )
    } else {
      data_plot <- ggplot2::ggplot() +
        ggplot2::geom_text(ggplot2::aes(x = 0.5, y = 0.5, label = "No data available"))
    }

    return(list(
      effect_2d = effect_plot,
      contour = contour_plot,
      data_distribution = data_plot,
      combined = gridExtra::grid.arrange(effect_plot, contour_plot, data_plot, ncol = 2)
    ))
  } else {
    # Variables not found in data
    return(list(
      effect_2d = ggplot2::ggplot() +
        ggplot2::geom_text(ggplot2::aes(
          x = 0.5, y = 0.5,
          label = paste(
            "Variables", paste(variables, collapse = ", "),
            "not found in data"
          )
        )),
      contour = ggplot2::ggplot(),
      data_distribution = ggplot2::ggplot()
    ))
  }
}

# CDI plot for parametric terms
cdi_plot_parametric <- function(term, variable, colour, smooth_info) {
  if (is.null(variable)) {
    variable <- smooth_info$variables[1]
  }

  # Use the original CDI approach for parametric terms
  return(private$cdi_plot_1d(term, variable, colour, 20, smooth_info))
}

# Prepare coefficient data for 1D terms
prepare_coefficient_data_1d <- function(term, variable, n_bins) {
  if (variable %in% names(self$data)) {
    var_data <- self$data[[variable]]

    if (is.numeric(var_data)) {
      # Create bins for continuous variables
      breaks <- pretty(var_data, n_bins)
      step <- breaks[2] - breaks[1]
      labels <- breaks[-length(breaks)] + step / 2
      breaks <- c(breaks, breaks[length(breaks)] + step)
      levels_var <- cut(var_data, breaks, labels = labels, include.lowest = TRUE)
    } else {
      levels_var <- var_data
    }

    # Create prediction data for each level
    unique_levels <- levels(levels_var)
    if (is.null(unique_levels)) unique_levels <- unique(levels_var)

    pred_data_list <- list()

    for (i in seq_along(unique_levels)) {
      level_val <- unique_levels[i]

      # Create representative data point for this level
      pred_point <- self$data[1, , drop = FALSE] # Start with first row

      if (is.numeric(var_data)) {
        pred_point[[variable]] <- as.numeric(as.character(level_val))
      } else {
        pred_point[[variable]] <- level_val
      }

      # Set other variables to their means/modes
      for (var in names(self$data)) {
        if (var != variable && var != self$response) {
          if (is.numeric(self$data[[var]])) {
            pred_point[[var]] <- mean(self$data[[var]], na.rm = TRUE)
          } else {
            pred_point[[var]] <- names(sort(table(self$data[[var]]), decreasing = TRUE))[1]
          }
        }
      }

      pred_data_list[[i]] <- pred_point
    }

    pred_data <- do.call(rbind, pred_data_list)

    # Get predictions
    tryCatch(
      {
        preds <- predict(self$model, newdata = pred_data, type = "terms", se.fit = TRUE)

        if (term %in% colnames(preds$fit)) {
          coeff_data <- data.frame(
            level = unique_levels,
            level_num = 1:length(unique_levels),
            coeff = preds$fit[, term],
            se = preds$se.fit[, term],
            stringsAsFactors = FALSE
          )
        } else {
          coeff_data <- data.frame(
            level = unique_levels,
            level_num = 1:length(unique_levels),
            coeff = rep(0, length(unique_levels)),
            se = rep(0, length(unique_levels)),
            stringsAsFactors = FALSE
          )
        }

        coeff_data$exp_coeff <- exp(coeff_data$coeff)
        coeff_data$exp_lower <- exp(coeff_data$coeff - coeff_data$se)
        coeff_data$exp_upper <- exp(coeff_data$coeff + coeff_data$se)

        return(coeff_data)
      },
      error = function(e) {
        # Return dummy data on error
        return(data.frame(
          level = unique_levels,
          level_num = 1:length(unique_levels),
          coeff = rep(0, length(unique_levels)),
          se = rep(0, length(unique_levels)),
          exp_coeff = rep(1, length(unique_levels)),
          exp_lower = rep(1, length(unique_levels)),
          exp_upper = rep(1, length(unique_levels)),
          stringsAsFactors = FALSE
        ))
      }
    )
  } else {
    # Variable not found, return empty data frame
    return(data.frame(
      level = character(0),
      level_num = numeric(0),
      coeff = numeric(0),
      se = numeric(0),
      exp_coeff = numeric(0),
      exp_lower = numeric(0),
      exp_upper = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
}

# Prepare distribution data for CDI plot
prepare_distribution_data <- function(term, variable, n_bins) {
  smooth_info <- private$get_smooth_info(self$focus)

  if (variable %in% names(self$data) &&
    smooth_info$type == "parametric" &&
    self$focus %in% names(self$data)) {
    var_data <- self$data[[variable]]

    if (is.numeric(var_data)) {
      breaks <- pretty(var_data, n_bins)
      step <- breaks[2] - breaks[1]
      labels <- breaks[-length(breaks)] + step / 2
      breaks <- c(breaks, breaks[length(breaks)] + step)
      levels_var <- cut(var_data, breaks, labels = labels, include.lowest = TRUE)
    } else {
      levels_var <- var_data
    }

    # Create distribution data
    distrib_data <- self$data %>%
      dplyr::mutate(term_level = levels_var) %>%
      dplyr::filter(!is.na(term_level)) %>%
      dplyr::group_by(term_level, !!rlang::sym(self$focus)) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(!!rlang::sym(self$focus)) %>%
      dplyr::mutate(
        total = sum(count),
        prop = count / total
      ) %>%
      dplyr::ungroup() %>%
      dplyr::rename(focus_level = !!rlang::sym(self$focus))

    return(distrib_data)
  } else {
    # Return empty data frame if variables not suitable
    return(data.frame(
      term_level = character(0),
      focus_level = character(0),
      count = numeric(0),
      total = numeric(0),
      prop = numeric(0)
    ))
  }
}
