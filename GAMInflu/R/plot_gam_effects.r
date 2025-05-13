#' Plot Model Effects for GAM
#'
#' Creates plots showing the effect of each predictor in a GAM model while holding
#' other variables at their median/mode/reference level. For random effects,
#' includes specialized visualizations showing the distribution and properties
#' of the random effect estimates.
#'
#' @param model A fitted GAM object (from mgcv::gam)
#' @param data The data frame used to fit the model
#' @param terms Character vector of terms to plot. If NULL (default), plots all terms
#' @param n_points Number of points to use for continuous variables
#' @param show_rug Logical. Whether to show rug plots for continuous variables
#' @param show_ci Logical. Whether to show confidence intervals
#' @param ci_level Confidence level (0-1)
#' @param plot_as_grid Logical. Whether to arrange plots in a grid (TRUE) or return a list of plots (FALSE)
#' @param smooth_method Method for smooth terms. Options: "predict" (uses predict.gam) or "plot_smooth" (uses plot.gam)
#' @param ncol Number of columns in the grid (if plot_as_grid=TRUE)
#' @param verbose Logical. Print progress messages?
#' @param re_style Style for random effects plots. Options: "panel" (default), "qqnorm", "caterpillar", "shrinkage", "combined"
#' @param interactive Logical. Whether to create interactive plots using plotly (if available)
#' @param predict_all Logical. Whether to predict for all terms even if not visualized
#'
#' @return A grid of plots (if plot_as_grid=TRUE) or a list of ggplot objects
#' @export
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' data(mtcars)
#' m <- gam(mpg ~ s(wt) + s(hp) + s(cyl, bs = "re"), data = mtcars)
#' plot_gam_effects(m, mtcars)
#' plot_gam_effects(m, mtcars, re_style = "qqnorm")
#' plot_gam_effects(m, mtcars, re_style = "caterpillar")
#' }
plot_gam_effects <- function(model, data, terms = NULL, n_points = 100,
                             show_rug = FALSE, show_ci = TRUE, ci_level = 0.95,
                             plot_as_grid = TRUE, smooth_method = "predict",
                             ncol = 2, verbose = TRUE, re_style = "panel",
                             interactive = FALSE, predict_all = TRUE) {
  # Check inputs
  if (!inherits(model, "gam")) {
    stop("Model must be a GAM object from mgcv package")
  }

  if (!is.data.frame(data)) {
    stop("Data must be a data frame")
  }
  re_style <- match.arg(re_style, c("panel", "qqnorm", "caterpillar", "shrinkage", "combined"))
  if (!is.logical(plot_as_grid)) {
    stop("plot_as_grid must be a logical value")
  }
  if (!is.logical(show_rug)) {
    stop("show_rug must be a logical value")
  }
  if (!is.logical(show_ci)) {
    stop("show_ci must be a logical value")
  }
  if (!is.numeric(n_points) || n_points <= 0) {
    stop("n_points must be a positive integer")
  }
  if (!is.numeric(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be between 0 and 1")
  }
  if (!is.numeric(ncol) || ncol <= 0) {
    stop("ncol must be a positive integer")
  }
  if (!is.logical(verbose)) {
    stop("verbose must be a logical value")
  }
  if (!is.logical(interactive)) {
    stop("interactive must be a logical value")
  }
  if (!is.logical(predict_all)) {
    stop("predict_all must be a logical value")
  }
  if (!is.character(smooth_method) || !smooth_method %in% c("predict", "plot_smooth")) {
    stop("smooth_method must be either 'predict' or 'plot_smooth'")
  }
  if (!is.character(re_style) || !re_style %in% c("panel", "qqnorm", "caterpillar", "shrinkage", "combined")) {
    stop("re_style must be one of 'panel', 'qqnorm', 'caterpillar', 'shrinkage', or 'combined'")
  }
  if (!is.null(terms) && !is.character(terms)) {
    stop("terms must be a character vector")
  }
  if (!is.null(terms) && length(terms) == 0) {
    stop("terms cannot be an empty vector")
  }
  if (!is.null(terms) && any(!terms %in% names(data))) {
    stop("Some terms are not present in the data")
  }
  if (!is.null(terms) && any(terms %in% names(data))) {
    terms <- terms[terms %in% names(data)]
  }
  if (length(terms) == 0) {
    stop("No valid terms to plot. Ensure terms match variable names in the model.")
  }
  if (length(terms) > 1 && plot_as_grid) {
    if (!is.numeric(ncol) || ncol <= 0) {
      stop("ncol must be a positive integer")
    }
  }
  if (length(terms) == 1 && plot_as_grid) {
    plot_as_grid <- FALSE
  }
  if (length(terms) > 1 && !plot_as_grid) {
    plot_as_grid <- TRUE
  }

  # Required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed to create plots")
  }

  # Extract response variable from model formula
  response_var <- as.character(stats::formula(model)[[2]])
  if (verbose) message("Response variable: ", response_var)

  # Extract all terms from model
  model_formula <- stats::formula(model)
  all_terms <- attr(stats::terms(model), "term.labels")

  # Extract variable names from smooth terms
  smooth_var_names <- sapply(model$smooth, function(smooth) smooth$term)
  smooth_var_names <- unique(unlist(smooth_var_names)) # Flatten list and remove duplicates

  # Identify random effect terms in the model
  random_terms_indices <- sapply(model$smooth, function(x) {
    inherits(x, "random.effect") || (is.character(x$bs) && x$bs == "re")
  })
  random_terms <- character(0)
  if (any(random_terms_indices)) {
    random_terms <- sapply(which(random_terms_indices), function(i) model$smooth[[i]]$term)
  }

  if (verbose && length(random_terms) > 0) {
    message("Identified random effect terms: ", paste(random_terms, collapse = ", "))
  }

  # Filter terms if specified
  if (!is.null(terms)) {
    # Match user-specified terms to variable names in the model
    matched_terms <- terms[terms %in% smooth_var_names]
    if (length(matched_terms) == 0) {
      stop("No valid terms to plot. Ensure terms match variable names in the model.")
    }
    terms <- matched_terms
  } else {
    terms <- smooth_var_names
  }

  if (verbose) {
    message("Creating effect plots for: ", paste(terms, collapse = ", "))
    if (length(random_terms) > 0) {
      message("Using random effect style: ", re_style)
    }
  }

  # Special case: if there are multiple random effects and we're using a combined view
  if (length(random_terms) > 1 && re_style == "combined" && plot_as_grid && all(random_terms %in% terms)) {
    # Create individual plots first, then add a combined view
    combined_re_viz <- TRUE
  } else {
    combined_re_viz <- FALSE
  }

  # Create prediction grid for each term
  plot_list <- list()

  for (term in terms) {
    if (verbose) message("Processing term: ", term)

    # Check if term is a smooth function
    is_smooth <- term %in% smooth_var_names

    # Check if term is a random effect
    is_random_effect <- term %in% random_terms

    # Find the smooth ID for this term
    smooth_id <- NULL
    if (is_smooth) {
      for (i in seq_along(model$smooth)) {
        if (term %in% model$smooth[[i]]$term) {
          smooth_id <- i
          break
        }
      }
    }

    # Extract variable name(s) from the term
    var_names <- term

    # Create prediction grid
    pred_grid <- create_prediction_grid(model, data, var_names, n_points, verbose)

    # Get predictions for this term
    if (is_smooth && smooth_method == "plot_smooth" && length(var_names) == 1) {
      # For univariate smooths, use plot.gam's built-in extraction
      plot_data <- extract_smooth_term(model, term, pred_grid, ci_level)
    } else {
      # For regular terms or complex smooths, use predict.gam
      plot_data <- predict_term_effect(model, term, pred_grid, ci_level, show_ci)
    }

    # Generate the appropriate plot based on term type
    if (is_random_effect) {
      # For random effects, create specialized visualizations
      p <- create_random_effect_plot(model, term, plot_data, data, var_names, response_var, smooth_id, show_rug, show_ci, re_style = re_style, verbose = verbose)
    } else if (is_smooth && length(var_names) >= 2) {
      # For spatial smooths (bivariate or higher)
      p <- create_spatial_plot(term, plot_data, data, var_names, show_rug, show_ci, color_palette = "viridis", contour_type = "raster")
    } else {
      # For regular terms
      p <- create_effect_plot(
        term, plot_data, data, var_names, is_smooth, show_rug, show_ci,
        response_var = response_var
      )
    }

    plot_list[[term]] <- p
  }

  # If we have multiple random effects and combined view is requested, add it
  if (combined_re_viz) {
    # Create combined random effects visualization
    p_combined <- create_combined_re_plot(model, random_terms, data, response_var)
    plot_list[["combined_random_effects"]] <- p_combined
  }

  # Convert to interactive plots if requested
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Package 'plotly' needed for interactive plots. Returning static plots instead.")
    } else {
      # Convert plots to plotly objects
      # Note: This won't work for grid objects from grid.arrange
      # So we skip those and add a warning
      for (i in seq_along(plot_list)) {
        if (inherits(plot_list[[i]], "ggplot")) {
          plot_list[[i]] <- plotly::ggplotly(plot_list[[i]])
        } else if (verbose) {
          message("  Note: Plot for term '", names(plot_list)[i], "' could not be converted to interactive format.")
        }
      }
    }
  }

  # Return result
  if (plot_as_grid && length(plot_list) > 1) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      warning("Package 'patchwork' needed for grid layout. Returning list of plots instead.")
      return(plot_list)
    }

    # Check if all plots are ggplot objects (or at least can be arranged by patchwork)
    can_arrange <- TRUE
    for (p in plot_list) {
      if (!inherits(p, "ggplot") && !inherits(p, "gtable")) {
        can_arrange <- FALSE
        break
      }
    }

    if (can_arrange) {
      # Arrange plots in a grid using patchwork
      result <- patchwork::wrap_plots(plot_list, ncol = ncol) + patchwork::plot_annotation()
      return(result)
    } else {
      warning("Some plots are not compatible with patchwork. Returning list of plots instead.")
      return(plot_list)
    }
  } else {
    return(plot_list)
  }
}

# Helper function for random effect diagnostics
compute_re_diagnostics <- function(model, smooth_id) {
  re_coefs <- coef(model)[model$smooth[[smooth_id]]$first.para:model$smooth[[smooth_id]]$last.para]
  diag_list <- list(
    n_groups = length(re_coefs),
    sd = stats::sd(re_coefs),
    mean = mean(re_coefs),
    min_max_ratio = ifelse(
      any(abs(re_coefs) > 0.0001),
      max(abs(re_coefs)) / min(abs(re_coefs[abs(re_coefs) > 0.0001])),
      NA
    ),
    shrinkage_factor = model$sp[smooth_id]
  )

  # Check normality with Shapiro-Wilk test (if we have enough data)
  if (length(re_coefs) >= 3 && length(re_coefs) <= 5000) {
    diag_list$shapiro_p <- stats::shapiro.test(re_coefs)$p.value
    diag_list$normality <- ifelse(diag_list$shapiro_p > 0.05, "Normal", "Non-normal")
  } else {
    diag_list$normality <- "Too many/few groups for test"
  }

  return(diag_list)
}

# Enhanced helper function to create plots with density panels for random effects
create_random_effect_plot <- function(model, term, plot_data, data, var_names,
                                      response_var, smooth_id, show_rug, show_ci,
                                      re_style = "panel", verbose = TRUE) {
  # Find the smooth corresponding to this term
  if (is.null(smooth_id)) {
    for (i in seq_along(model$smooth)) {
      if (model$smooth[[i]]$label == term || grepl(paste0("*\\(", term, "\\)*"), model$smooth[[i]]$label)) {
        smooth_id <- i
        break
      }
    }
  }

  if (is.null(smooth_id)) {
    warning("Could not find smooth term: ", term, ". Falling back to standard plot.")
    return(create_effect_plot(term, plot_data, data, var_names, TRUE, show_rug, show_ci))
  }

  # Extract the random effect coefficients
  smooth <- model$smooth[[smooth_id]]
  coef_indices <- smooth$first.para:smooth$last.para
  random_coefs <- coef(model)[coef_indices]

  # Compute diagnostics
  diagnostics <- compute_re_diagnostics(model, smooth_id)

  # For multi-variable terms, only plot first variable
  x_var <- var_names[1]

  # Create the main effect plot
  if (is.numeric(plot_data[[x_var]])) {
    # Continuous variable
    p_main <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[["y"]])) +
      ggplot2::geom_line(aes(colour = "blue", alpha = 0.5)) +
      ggplot2::labs(x = x_var, y = "Partial effects")

    if (show_ci && "lower" %in% names(plot_data) && "upper" %in% names(plot_data)) {
      p_main <- p_main + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]], colour = "blue"), alpha = 0.2)
    }

    if (show_rug && x_var %in% names(data)) {
      p_main <- p_main + ggplot2::geom_rug(data = data, ggplot2::aes(x = .data[[x_var]], y = NULL), sides = "b", alpha = 0.2)
    }
  } else {
    # Categorical variable
    p_main <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[["y"]])) +
      ggplot2::geom_point(aes(colour = "blue", alpha = 0.5)) +
      ggplot2::labs(x = x_var, y = "Partial effects")

    if (show_ci && "lower" %in% names(plot_data) && "upper" %in% names(plot_data)) {
      p_main <- p_main + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]), width = 0.2)
    }
    # For factors, adjust x-axis
    if (is.factor(plot_data[[x_var]])) {
      p_main <- p_main + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
  }

  # Create diagnostics annotation
  diag_text <- paste(
    "Random Effect Diagnostics:",
    paste("N groups:", diagnostics$n_groups),
    paste("SD:", round(diagnostics$sd, 4)),
    paste("Normality:", diagnostics$normality),
    paste("Shrinkage factor:", round(diagnostics$shrinkage_factor, 4)),
    sep = "\n"
  )

  # Create appropriate secondary plot based on re_style
  if (re_style == "panel") {
    # Density panel
    density_df <- data.frame(effect = random_coefs)

    p_secondary <- ggplot2::ggplot(density_df, ggplot2::aes(x = effect)) +
      ggplot2::geom_density(fill = "lightblue", alpha = 0.7) +
      ggplot2::geom_rug(alpha = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::labs(x = "Effect", y = "Density") +
      ggplot2::coord_flip()
  } else if (re_style == "qqnorm") {
    if (verbose) message("  Creating QQ-plot for random effects")

    # Ensure random_coefs is a numeric vector
    density_df <- data.frame(effect = random_coefs)

    # Create the Q-Q plot
    p_secondary <- ggplot2::ggplot(density_df, ggplot2::aes(sample = effect)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line() +
      ggplot2::labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
  } else if (re_style == "caterpillar") {
    if (verbose) message("  Creating caterpillar plot for random effects")

    # Calculate standard errors for random effects
    se_val <- sqrt(1 / model$sp[smooth_id]) # Approximate SE based on smoothing parameter

    # Prepare data for caterpillar plot
    caterpillar_data <- data.frame(
      level = factor(seq_along(random_coefs)),
      effect = random_coefs,
      se = se_val
    )

    # Sort by effect size
    caterpillar_data <- caterpillar_data[order(caterpillar_data$effect), ]
    caterpillar_data$level <- factor(caterpillar_data$level, levels = caterpillar_data$level)

    # Add confidence intervals
    caterpillar_data$lower <- caterpillar_data$effect - 1.96 * caterpillar_data$se
    caterpillar_data$upper <- caterpillar_data$effect + 1.96 * caterpillar_data$se

    # Create the caterpillar plot
    p_secondary <- ggplot2::ggplot(caterpillar_data, ggplot2::aes(x = level, y = effect)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2, alpha = 0.5) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(x = "Groups (sorted by effect)", y = "Random Effect") +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  } else if (re_style == "shrinkage") {
    # Need to calculate raw group means for comparison
    # First identify which variable in the data corresponds to the grouping factor
    group_var <- var_names[1]

    # Check if we can extract group means
    if (group_var %in% names(data) && response_var %in% names(data)) {
      # Calculate group means
      if (requireNamespace("dplyr", quietly = TRUE)) {
        if (verbose) message("  Creating shrinkage plot for random effects")

        group_data <- data %>%
          dplyr::group_by(!!dplyr::sym(group_var)) %>%
          dplyr::summarize(raw_mean = mean(!!dplyr::sym(response_var), na.rm = TRUE))

        # Match with random effects
        if (nrow(group_data) == length(random_coefs)) {
          shrinkage_data <- data.frame(
            raw_mean = group_data$raw_mean,
            re_estimate = random_coefs
          )

          # Create shrinkage plot
          p_secondary <- ggplot2::ggplot(shrinkage_data, ggplot2::aes(x = raw_mean, y = re_estimate)) +
            ggplot2::geom_point() +
            ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
            ggplot2::labs(x = "Raw Group Means", y = "Random Effect Estimates")
        } else {
          # Fallback to panel if mismatch
          if (verbose) {
            warning("Cannot create shrinkage plot: mismatch between group levels and random effects.")
            message("  Falling back to panel style")
          }
          re_style <- "panel"
        }
      } else {
        # Fallback if dplyr not available
        if (verbose) {
          warning("Package 'dplyr' needed for shrinkage plots. Falling back to panel style.")
        }
        re_style <- "panel"
      }
    } else {
      # Fallback if required variables not available
      if (verbose) {
        warning("Cannot create shrinkage plot: required variables not in data.")
        message("  Falling back to panel style")
      }
      re_style <- "panel"
    }
  } else if (re_style == "panel") {
    density_df <- data.frame(effect = random_coefs)

    p_secondary <- ggplot2::ggplot(density_df, ggplot2::aes(x = effect)) +
      ggplot2::geom_density(fill = "lightblue", alpha = 0.7) +
      ggplot2::geom_rug(alpha = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::labs(x = "Effect", y = "Density") +
      ggplot2::coord_flip()
  } else {
    # Default to density panel if unknown style
    if (verbose) message("  Unknown re_style: ", re_style, ". Using panel style.")
    density_df <- data.frame(effect = random_coefs)

    p_secondary <- ggplot2::ggplot(density_df, ggplot2::aes(x = effect)) +
      ggplot2::geom_density(fill = "lightblue", alpha = 0.7) +
      ggplot2::geom_rug(alpha = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::labs(x = "Effect", y = "Density", ) +
      ggplot2::coord_flip()
  }

  # Combine the main and secondary plots
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' needed for combined plots. Returning main plot only.")
    return(p_main)
  }

  combined_plot <- p_main + p_secondary +
    patchwork::plot_layout(ncol = 2, widths = c(3, 2))
  return(combined_plot)
}

# Function to create a combined visualization of all random effects
create_combined_re_plot <- function(model, random_terms, data, response_var) {
  # Initialize empty data frame for all random effects
  all_re_data <- data.frame()

  # Extract random effects for each term
  for (term in random_terms) {
    # Find the smooth ID
    smooth_id <- NULL
    for (i in seq_along(model$smooth)) {
      if (model$smooth[[i]]$label == term) {
        smooth_id <- i
        break
      }
    }

    if (!is.null(smooth_id)) {
      # Extract coefficients
      coef_indices <- model$smooth[[smooth_id]]$first.para:model$smooth[[smooth_id]]$last.para
      random_coefs <- coef(model)[coef_indices]

      # Add to data frame
      term_data <- data.frame(
        effect = random_coefs,
        term = term
      )

      all_re_data <- rbind(all_re_data, term_data)
    }
  }

  # Create density plot for all random effects
  if (nrow(all_re_data) > 0) {
    p_combined <- ggplot2::ggplot(all_re_data, ggplot2::aes(x = effect, fill = term)) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::facet_wrap(~term, scales = "free") +
      ggplot2::labs(x = "Partial effects", y = "Density")
    print(p_combined)
    return(p_combined)
  } else {
    # Return NULL if no data available
    return(NULL)
  }
}

# Helper function to create prediction grid
create_prediction_grid <- function(model, data, var_names, n_points, verbose) {
  # Start with all variables at reference/median values
  pred_grid <- data.frame(row.names = 1:n_points)

  # Get all variables from the model
  model_vars <- all.vars(stats::formula(model)[-2]) # Exclude response

  # Add reference values for all model variables
  for (var in model_vars) {
    if (var %in% names(data)) {
      if (is.factor(data[[var]]) || is.character(data[[var]])) {
        # For factors, use the reference level or most common level
        if (is.factor(data[[var]])) {
          ref_level <- levels(data[[var]])[1] # Reference level
        } else {
          # For character, convert to factor first
          tab <- table(data[[var]])
          ref_level <- names(tab)[which.max(tab)] # Most common value
        }
        pred_grid[[var]] <- rep(ref_level, n_points)
      } else if (is.numeric(data[[var]])) {
        # For numeric, use the median
        pred_grid[[var]] <- rep(stats::median(data[[var]], na.rm = TRUE), n_points)
      } else {
        # For other types, use the first value
        pred_grid[[var]] <- rep(data[[var]][1], n_points)
      }
    }
  }

  # Now set the sequence for variables of interest
  for (var in var_names) {
    if (var %in% names(data)) {
      if (is.factor(data[[var]]) || is.character(data[[var]])) {
        # For factors, use all levels
        if (is.factor(data[[var]])) {
          levels_to_use <- intersect(levels(data[[var]]), levels(data[[var]]))
          pred_grid[[var]] <- factor(rep(levels_to_use, length.out = n_points), levels = levels(data[[var]]))
        } else {
          # For character variables, use unique values
          unique_vals <- unique(data[[var]])
          pred_grid[[var]] <- rep(unique_vals, length.out = n_points)
        }
      } else if (is.numeric(data[[var]])) {
        # For numeric, use sequence over range
        x_range <- range(data[[var]], na.rm = TRUE)
        pred_grid[[var]] <- seq(from = x_range[1], to = x_range[2], length.out = n_points)
      } else {
        if (verbose) message("  Warning: Variable ", var, " has unsupported type. Using first value.")
      }
    } else {
      if (verbose) message("  Warning: Variable ", var, " not found in data.")
    }
  }

  return(pred_grid)
}

# Helper function to extract smooth term using plot.gam
extract_smooth_term <- function(model, term, pred_grid, ci_level) {
  # Extract the variable name from the smooth term
  var_name <- gsub("^s\\(|^te\\(|^ti\\(|\\)$|,.+\\)$", "", term)

  # Get smooths list from model
  smooths <- model$smooth

  # Find the appropriate smooth
  smooth_id <- NULL
  for (i in seq_along(smooths)) {
    if (smooths[[i]]$label == term || # Check label
      # Or check if the term name is in the smooth
      grepl(paste0("^", term, "$"), smooths[[i]]$label)) {
      smooth_id <- i
      break
    }
  }

  if (is.null(smooth_id)) {
    stop("Could not find smooth term: ", term)
  }

  # Extract the smooth
  smooth <- smooths[[smooth_id]]

  # Use mgcv:::plot.gam to extract the partial residuals
  tmp <- tempfile()
  pdf(tmp)
  # Handle potential errors
  p_obj <- tryCatch(
    {
      p_obj <- plot(model, select = smooth_id, se = (ci_level > 0))
      p_obj
    },
    error = function(e) {
      warning("Error in plot.gam for term ", term, ": ", e$message)
      NULL
    }
  )
  unlink(tmp)

  # If plot.gam fails, fall back to predict method
  if (is.null(p_obj)) {
    return(predict_term_effect(model, term, pred_grid, ci_level, TRUE))
  }

  # Extract the data from the plot object
  x_var <- p_obj[[1]]$x
  y_var <- p_obj[[1]]$fit

  if (ci_level > 0) {
    # Confidence interval
    ci_mult <- stats::qnorm((1 + ci_level) / 2)
    y_lower <- y_var - ci_mult * p_obj[[1]]$se
    y_upper <- y_var + ci_mult * p_obj[[1]]$se
  } else {
    y_lower <- y_upper <- NULL
  }

  # Create a data frame with the results
  plot_data <- data.frame(
    x = x_var,
    y = y_var,
    lower = y_lower,
    upper = y_upper
  )
  names(plot_data)[1] <- var_name

  return(plot_data)
}

# Helper function to predict term effect using predict.gam
predict_term_effect <- function(model, term, pred_grid, ci_level, show_ci) {
  # Predict for all data points
  if (show_ci) {
    preds <- stats::predict(model, newdata = pred_grid, type = "link", se.fit = TRUE)
    fit <- preds$fit
    se <- preds$se.fit

    # Confidence interval
    ci_mult <- stats::qnorm((1 + ci_level) / 2)
    lower <- fit - ci_mult * se
    upper <- fit + ci_mult * se

    # Transform to response scale if needed
    if (inherits(model$family, "family")) {
      fit <- model$family$linkinv(fit)
      lower <- model$family$linkinv(lower)
      upper <- model$family$linkinv(upper)
    }

    preds_df <- data.frame(pred_grid, y = fit, lower = lower, upper = upper)
  } else {
    fit <- stats::predict(model, newdata = pred_grid, type = "response")
    preds_df <- data.frame(pred_grid, y = fit)
  }

  return(preds_df)
}

# Helper function to create the plot
create_effect_plot <- function(term, plot_data, data, var_names, is_smooth, show_rug, show_ci,
                               response_var = NULL) {
  # For multi-variable terms, only plot first variable
  x_var <- var_names[1]

  # Set y-axis label
  y_lab <- "Partial effects"

  # Handle different types of variables
  if (is.numeric(plot_data[[x_var]])) {
    # Continuous variable
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[["y"]])) +
      ggplot2::geom_line(colour = "blue", alpha = 0.8) +
      ggplot2::labs(x = x_var, y = y_lab)

    if (show_ci && "lower" %in% names(plot_data) && "upper" %in% names(plot_data)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
        fill = "blue", alpha = 0.3
      )
    }

    if (show_rug && x_var %in% names(data)) {
      p <- p + ggplot2::geom_rug(
        data = data,
        ggplot2::aes(x = .data[[x_var]], y = NULL),
        sides = "b",
        alpha = 0.2
      )
    }
  } else {
    # Categorical variable
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(
      x = .data[[x_var]],
      y = .data[["y"]],
      fill = .data[[group_var]],
      colour = .data[[group_var]]
    )) +
      ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.4) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.9), alpha = 0.4) +
      ggplot2::labs(x = x_var, y = y_lab, fill = group_var, colour = group_var)

    if (show_ci && "lower" %in% names(plot_data) && "upper" %in% names(plot_data)) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
        position = ggplot2::position_dodge(width = 0.9),
        width = 0.2,
        colour = "blue",
        alpha = 0.3
      )
    }

    # For factors, adjust x-axis
    if (is.factor(plot_data[[x_var]])) {
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }

    if (color_palette == "viridis" && requireNamespace("viridisLite", quietly = TRUE)) {
      p <- p + ggplot2::scale_fill_viridis_d() + ggplot2::scale_colour_viridis_d()
    }
  }

  return(p)
}

# Helper function for spatial smooths - can be used for MRF random effects too
create_spatial_plot <- function(term, plot_data, data, var_names, show_rug, show_ci, color_palette, contour_type) {
  # Need at least 2 variables for spatial plot
  if (length(var_names) < 2) {
    stop("Spatial plot requires at least two variables.")
  }

  var1 <- var_names[1]
  var2 <- var_names[2]

  # Check if variables are numeric
  if (!is.numeric(plot_data[[var1]]) || !is.numeric(plot_data[[var2]])) {
    # Fall back to interaction plot for non-numeric variables
    warning("Non-numeric variables in spatial term. Falling back to standard plot.")
    p <- create_effect_plot(term, plot_data, data, var_names, TRUE, show_rug, show_ci)
    return(p)
  }

  # Create the base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[var1]], y = .data[[var2]]))

  # Add contours based on type
  if (contour_type == "raster") {
    # Filled contour with color scale
    p <- p + ggplot2::geom_raster(ggplot2::aes(fill = .data[["y"]])) +
      ggplot2::geom_contour(ggplot2::aes(z = .data[["y"]]), color = "white", alpha = 0.5, linewidth = 0.2)

    # Apply color scale
    if (color_palette == "viridis" && requireNamespace("viridisLite", quietly = TRUE)) {
      p <- p + ggplot2::scale_fill_viridis_c(option = "plasma")
    } else if (requireNamespace("RColorBrewer", quietly = TRUE) &&
      color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      # Use RColorBrewer sequential palette
      p <- p + ggplot2::scale_fill_distiller(palette = color_palette, direction = 1)
    }
  } else {
    # Line contours
    p <- p + ggplot2::geom_contour(ggplot2::aes(z = .data[["y"]]), color = "darkblue") +
      ggplot2::scale_color_viridis_c(option = "plasma")
  }

  # Add rug plots if requested
  if (show_rug) {
    p <- p + ggplot2::geom_rug(
      data = data,
      ggplot2::aes(x = .data[[var1]], y = .data[[var2]]),
      alpha = 0.1,
      size = 0.1
    )
  }

  # Add labels and styling
  p <- p + ggplot2::labs(x = var1, y = var2, fill = if (contour_type == "raster") "Effect" else NULL, color = if (contour_type != "raster") "Effect" else NULL)

  return(p)
}

# Helper function for handling interaction terms
create_interaction_plot <- function(term, plot_data, data, var_names, is_smooth,
                                    show_rug, show_ci, color_palette = "viridis") {
  # Needs at least 2 variables
  if (length(var_names) < 2) {
    stop("Interaction plot requires at least two variables.")
  }

  # For interaction plots, we use the first two variables
  x_var <- var_names[1]
  group_var <- var_names[2]

  # Create the appropriate plot based on variable types
  if (is.numeric(plot_data[[x_var]])) {
    # Continuous x variable
    if (is.factor(plot_data[[group_var]]) || is.character(plot_data[[group_var]])) {
      if (!is.factor(plot_data[[group_var]])) {
        plot_data[[group_var]] <- factor(plot_data[[group_var]])
      }

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[["y"]], color = .data[[group_var]], group = .data[[group_var]])) +
        ggplot2::geom_line(alpha = 0.8) +
        ggplot2::labs(x = x_var, y = "Partial Effect", color = group_var)

      if (show_ci && "lower" %in% names(plot_data) && "upper" %in% names(plot_data)) {
        p <- p + ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]], fill = .data[[group_var]]),
          alpha = 0.3, color = NA
        )
      }

      if (show_rug && x_var %in% names(data)) {
        p <- p + ggplot2::geom_rug(data = data, ggplot2::aes(x = .data[[x_var]], y = NULL, color = .data[[group_var]]), sides = "b", alpha = 0.2)
      }

      if (color_palette == "viridis" && requireNamespace("viridisLite", quietly = TRUE)) {
        p <- p +
          ggplot2::scale_color_viridis_d() +
          ggplot2::scale_fill_viridis_d(alpha = 0.1)
      }
    } else if (is.numeric(plot_data[[group_var]])) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[[group_var]])) +
        ggplot2::geom_raster(ggplot2::aes(fill = .data[["y"]])) +
        ggplot2::geom_contour(ggplot2::aes(z = .data[["y"]]), color = "white", alpha = 0.5) +
        ggplot2::labs(x = x_var, y = group_var, fill = "Effect")

      if (color_palette == "viridis" && requireNamespace("viridisLite", quietly = TRUE)) {
        p <- p + ggplot2::scale_fill_viridis_c(option = "plasma")
      }
    }
  } else {
    # Categorical x variable
    if (is.factor(plot_data[[group_var]]) || is.character(plot_data[[group_var]])) {
      if (!is.factor(plot_data[[x_var]])) {
        plot_data[[x_var]] <- factor(plot_data[[x_var]])
      }
      if (!is.factor(plot_data[[group_var]])) {
        plot_data[[group_var]] <- factor(plot_data[[group_var]])
      }

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(
        x = .data[[x_var]],
        y = .data[["y"]],
        fill = .data[[group_var]],
        colour = .data[[group_var]]
      )) +
        ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.4) +
        ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.9), alpha = 0.4) +
        ggplot2::labs(x = x_var, y = "Partial Effect", fill = group_var, colour = group_var)

      if (show_ci && "lower" %in% names(plot_data) && "upper" %in% names(plot_data)) {
        p <- p + ggplot2::geom_errorbar(
          ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
          position = ggplot2::position_dodge(width = 0.9),
          width = 0.2,
          colour = "blue",
          alpha = 0.3
        )
      }

      if (color_palette == "viridis" && requireNamespace("viridisLite", quietly = TRUE)) {
        p <- p + ggplot2::scale_fill_viridis_d() + ggplot2::scale_colour_viridis_d()
      }
    } else if (is.numeric(plot_data[[group_var]])) {
      if (!is.factor(plot_data[[x_var]])) {
        plot_data[[x_var]] <- factor(plot_data[[x_var]])
      }

      breaks <- pretty(range(plot_data[[group_var]]), n = 5)
      plot_data$group_binned <- cut(plot_data[[group_var]], breaks = breaks)

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[["y"]], fill = .data$group_binned)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge", colour = "blue", alpha = 0.4) +
        ggplot2::labs(x = x_var, y = "Partial Effect", fill = group_var)

      if (color_palette == "viridis" && requireNamespace("viridisLite", quietly = TRUE)) {
        p <- p + ggplot2::scale_fill_viridis_d()
      }
    }
  }

  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(p)
}

# Function to extract a specific smooth term's contribution
extract_smooth_contribution <- function(model, smooth_id, data = NULL) {
  if (is.null(data)) {
    data <- model$model
  }

  # Get the smooth
  smooth <- model$smooth[[smooth_id]]

  # Create a prediction with only this smooth
  exclude <- rep(TRUE, length(coef(model)))
  exclude[smooth$first.para:smooth$last.para] <- FALSE

  # Create a zero model first
  zero_predictions <- predict(model,
    newdata = data, type = "link",
    exclude = seq_along(coef(model))
  )

  # Now predict with just this smooth
  smooth_predictions <- predict(model,
    newdata = data, type = "link",
    exclude = which(exclude)
  )

  # The contribution is the difference
  contribution <- smooth_predictions - zero_predictions

  return(contribution)
}
