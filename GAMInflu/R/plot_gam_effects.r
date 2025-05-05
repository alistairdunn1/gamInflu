#' Plot Model Effects for GAM
#'
#' Creates plots showing the effect of each predictor in a GAM model while holding
#' other variables at their median/mode/reference level.
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
#'
#' @return A grid of plots (if plot_as_grid=TRUE) or a list of ggplot objects
#' @export
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' data(mtcars)
#' m <- gam(mpg ~ s(wt) + s(hp) + factor(cyl), data = mtcars)
#' plot_gam_effects(m, mtcars)
#' }
plot_gam_effects <- function(model, data, terms = NULL, n_points = 100, 
                             show_rug = TRUE, show_ci = TRUE, ci_level = 0.95,
                             plot_as_grid = TRUE, smooth_method = "predict",
                             ncol = 2, verbose = TRUE) {
  
  # Check inputs
  if (!inherits(model, "gam")) {
    stop("Model must be a GAM object from mgcv package")
  }
  
  if (!is.data.frame(data)) {
    stop("Data must be a data frame")
  }
  
  # Required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed to create plots")
  }
  
  # Extract all terms from model
  model_formula <- stats::formula(model)
  all_terms <- attr(stats::terms(model), "term.labels")
  
  # Filter terms if specified
  if (!is.null(terms)) {
    if (!all(terms %in% all_terms)) {
      invalid_terms <- setdiff(terms, all_terms)
      warning("Some specified terms not found in model: ", 
              paste(invalid_terms, collapse = ", "))
      terms <- intersect(terms, all_terms)
    }
    if (length(terms) == 0) {
      stop("No valid terms to plot")
    }
  } else {
    terms <- all_terms
  }
  
  if (verbose) {
    message("Creating effect plots for: ", paste(terms, collapse = ", "))
  }
  
  # Create prediction grid for each term
  plot_list <- list()
  
  for (term in terms) {
    if (verbose) message("Processing term: ", term)
    
    # Check if term is a smooth function
    is_smooth <- grepl("^s\\(|^te\\(|^ti\\(", term)
    
    # Extract variable name(s) from the term
    if (is_smooth) {
      # Extract variable name from smooth term - simple extraction, might need improvement for complex cases
      var_name <- gsub("^s\\(|^te\\(|^ti\\(|\\)$|,.+\\)$", "", term)
      # Handle potential by= variable
      if (grepl("by\\s*=", term)) {
        by_var <- gsub(".*by\\s*=\\s*([^,)]+).*", "\\1", term)
        var_name <- c(var_name, by_var)
        if (verbose) message("  Found smooth with 'by' variable: ", by_var)
      }
      var_names <- unlist(strsplit(var_name, "\\s*,\\s*"))
    } else {
      # For regular terms including interactions (a:b)
      var_names <- unlist(strsplit(term, "\\s*:\\s*"))
    }
    
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
    
    # Generate plot
    p <- create_effect_plot(term, plot_data, data, var_names, is_smooth, show_rug, show_ci)
    plot_list[[term]] <- p
  }
  
  # Return result
  if (plot_as_grid && length(plot_list) > 1) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      warning("Package 'patchwork' needed for grid layout. Returning list of plots instead.")
      return(plot_list)
    }
    # Arrange plots in a grid
    result <- patchwork::wrap_plots(plot_list, ncol = ncol) + 
      patchwork::plot_annotation(
        title = "Model Effects Plot",
        subtitle = paste("Response variable:", deparse(model_formula[[2]]))
      )
    return(result)
  } else {
    return(plot_list)
  }
}

# Helper function to create prediction grid
create_prediction_grid <- function(model, data, var_names, n_points, verbose) {
  # Start with all variables at reference/median values
  pred_grid <- data.frame(row.names = 1:n_points)
  
  # Get all variables from the model
  model_vars <- all.vars(stats::formula(model)[-2])  # Exclude response
  
  # Add reference values for all model variables
  for (var in model_vars) {
    if (var %in% names(data)) {
      if (is.factor(data[[var]]) || is.character(data[[var]])) {
        # For factors, use the reference level or most common level
        if (is.factor(data[[var]])) {
          ref_level <- levels(data[[var]])[1]  # Reference level
        } else {
          # For character, convert to factor first
          tab <- table(data[[var]])
          ref_level <- names(tab)[which.max(tab)]  # Most common value
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
          levels_to_use <- levels(data[[var]])
          pred_grid[[var]] <- factor(rep(levels_to_use, length.out = n_points), 
                                     levels = levels(data[[var]]))
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
  p_obj <- tryCatch({
    p_obj <- plot(model, select = smooth_id, se = (ci_level > 0))
    p_obj
  }, error = function(e) {
    warning("Error in plot.gam for term ", term, ": ", e$message)
    NULL
  })
  dev.off()
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
create_effect_plot <- function(term, plot_data, data, var_names, is_smooth, show_rug, show_ci) {
  # For multi-variable terms, only plot first variable
  x_var <- var_names[1]
  
  # Handle different types of variables
  if (is.numeric(plot_data[[x_var]])) {
    # Continuous variable
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[["y"]])) +
      ggplot2::geom_line() +
      ggplot2::labs(
        title = paste("Effect of", term),
        x = x_var,
        y = "Partial Effect"
      )
    
    if (show_ci && "lower" %in% names(plot_data) && "upper" %in% names(plot_data)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
        alpha = 0.2
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
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[["y"]])) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(
        title = paste("Effect of", term),
        x = x_var,
        y = "Partial Effect"
      )
    
    if (show_ci && "lower" %in% names(plot_data) && "upper" %in% names(plot_data)) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
        width = 0.2
      )
    }
    
    # For factors, adjust x-axis
    if (is.factor(plot_data[[x_var]])) {
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
  }
  
  # Additional styling
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}

