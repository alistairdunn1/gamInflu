#' A coefficient-distribution-influence (CDI) plot for a model term
#'
#' @param obj The influence object
#' @param term The term to plot
#' @param variable Optional variable to use for the term
#' @param ... Additional arguments
#' @return A combined plot
#' @export
cdi_plot <- function(obj, term, variable = NULL, ...) {
  if (is.null(obj$preds) || is.null(obj$influences)) {
    stop("Call calc() before plotting")
  }
  
  # Check if this is an interaction term
  is_interaction <- term %in% names(obj$interactions)
  
  # For interaction terms, we need special handling
  if (is_interaction) {
    interaction_info <- obj$interactions[[term]]
    
    # If no variable is specified, use the first non-focus variable in the interaction
    if (is.null(variable)) {
      for (var in interaction_info$vars) {
        if (var != obj$focus) {
          variable <- var
          break
        }
      }
      
      # If all variables are the focus, just use the first variable
      if (is.null(variable)) {
        variable <- interaction_info$vars[1]
      }
    }
  } else {
    # For regular terms, find the variable that matches the term
    if (is.null(variable)) {
      for (name in names(obj$preds)) {
        # Look for the variable in the term
        pattern <- paste0("([^[:alnum:]_])+", name, "([^[:alnum:]_])+")
        term_with_brackets <- paste0("(", term, ")")
        match <- grep(pattern, term_with_brackets)
        if (length(match) > 0) {
          variable <- name
          break
        }
      }
      
      # If still not found, try exact match
      if (is.null(variable) && term %in% names(obj$preds)) {
        variable <- term
      }
    }
  }
  
  if (is.null(variable)) {
    # For smooth terms, use the term itself as variable
    if (grepl("^s\\(", term) || grepl("^ti\\(", term)) {
      # Extract variable name from smooth term 
      var_match <- regmatches(term, regexec("\\(([^,)]+)", term))
      if (length(var_match) > 0 && length(var_match[[1]]) > 1) {
        variable <- var_match[[1]][2]
      } else {
        # Fallback to term itself
        variable <- term
      }
    } else {
      stop('Unable to find a matching variable for term "', term, '". Please specify using the variable argument.')
    }
  }
  
  # Extract levels or create prediction grid
  if (is_interaction) {
    # For interactions, create a prediction grid
    pred_data <- create_prediction_grid(obj, term)
    levels <- pred_data[[variable]]
  } else {
    # For regular terms
    levels <- obj$preds[[variable]]
  }
  
  # Numeric terms are cut into factors
  if (is.numeric(levels)) {
    breaks <- pretty(levels, 30)
    step <- breaks[2] - breaks[1]
    labels <- breaks + step / 2
    breaks <- c(breaks, breaks[length(breaks)] + step)
    levels <- cut(levels, breaks, labels = labels, include.lowest = TRUE)
  }
  
  # Reorder levels according to coefficients if necessary
  if (term %in% names(obj$orders) && obj$orders[[term]] == "coef") {
    fit_term <- paste0("fit.", term)
    if (fit_term %in% names(obj$preds)) {
      coeffs <- aggregate(obj$preds[[fit_term]], list(levels), mean)
      names(coeffs) <- c("term", "coeff")
      coeffs <- coeffs[order(coeffs$coeff), ]
      levels <- factor(levels, levels = coeffs$term, ordered = TRUE)
    }
  }
  
  # Get term coefficients and standard errors
  # For interactions or complex terms, we may need to use predictions
  if (is_interaction || !paste0("fit.", term) %in% names(obj$preds)) {
    # Create a prediction grid for this term
    pred_data <- create_prediction_grid(obj, term)
    
    # Get predictions
    preds <- predict(obj$model, newdata = pred_data, type = "terms", terms = term, se.fit = TRUE)
    
    # Extract coefficients
    coeffs_df <- data.frame(
      term = pred_data[[variable]],
      coeff = preds$fit,
      se = preds$se.fit
    )
    
    # Aggregate by variable
    coeffs <- aggregate(coeffs_df[, c("coeff", "se")], list(term = coeffs_df$term), mean)
  } else {
    # For regular terms, use standard approach
    fit_term <- paste0("fit.", term)
    se_term <- paste0("se.fit.", term)
    
    coeffs <- aggregate(obj$preds[, c(fit_term, se_term)], list(term = levels), mean)
    names(coeffs) <- c("term", "coeff", "se")
  }
  
  # Prepare coefficient data for plotting
  coeffs <- coeffs %>%
    mutate(
      lower = coeff - se,
      upper = coeff + se,
      term_int = as.integer(term)
    )
  
  # Distribution plot data - handle interactions
  focus_var <- obj$preds[[obj$focus]]
  
  if (is_interaction) {
    # For interactions, use the prediction grid
    pred_data <- create_prediction_grid(obj, term)
    
    # Count occurrences of each combination
    distrs <- as.data.frame(table(pred_data[[variable]], pred_data[[obj$focus]]))
    names(distrs) <- c("term", "focus", "count")
    
    # Calculate proportions
    distrs <- merge(distrs, aggregate(list(total = distrs$count), list(focus = distrs$focus), sum))
    distrs <- distrs %>%
      mutate(
        prop = count / total,
        term_int = as.integer(term),
        focus_int = as.integer(focus)
      )
  } else {
    # Standard approach for regular terms
    distrs <- aggregate(rep(1, length(focus_var)), list(term = levels, focus = focus_var), length)
    names(distrs) <- c("term", "focus", "count")
    distrs <- merge(distrs, aggregate(list(total = distrs$count), list(focus = distrs$focus), sum))
    distrs <- distrs %>%
      mutate(
        prop = count / total,
        term_int = as.integer(term),
        focus_int = as.integer(focus)
      )
  }
  
  # Influence plot data - handle interactions
  if (is_interaction) {
    # For interaction terms, use the influences calculated in calc()
    if (term %in% names(obj$influences)) {
      influence_data <- obj$influences %>%
        dplyr::select(level, all_of(term)) %>%
        mutate(
          influence = exp(get(term)), 
          level_int = as.integer(level)
        )
    } else {
      # If the influence isn't directly available, calculate it
      pred_data <- create_prediction_grid(obj, term)
      preds <- predict(obj$model, newdata = pred_data, type = "terms", terms = term)
      
      # Aggregate by focus level
      infl <- aggregate(
        list(value = preds),
        by = list(level = pred_data[[obj$focus]]),
        mean
      )
      
      influence_data <- infl %>%
        mutate(
          influence = exp(value),
          level_int = as.integer(level)
        )
    }
  } else {
    # Standard approach for regular terms
    if (term %in% names(obj$influences)) {
      influence_data <- obj$influences %>%
        dplyr::select(level, all_of(term)) %>%
        mutate(
          influence = exp(get(term)), 
          level_int = as.integer(level)
        )
    } else {
      # Create placeholder influence data
      influence_data <- data.frame(
        level = levels(obj$data[[obj$focus]]),
        influence = 1,
        level_int = as.integer(levels(obj$data[[obj$focus]]))
      )
    }
  }
  
  # Create a nicer title based on term type
  if (is_interaction) {
    interaction_vars <- obj$interactions[[term]]$vars
    if (obj$interactions[[term]]$type %in% c("tensor.smooth", "tp", "ti")) {
      plot_title <- paste0("CDI Plot for ti(", paste(interaction_vars, collapse=","), ")")
    } else {
      plot_title <- paste0("CDI Plot for ", paste(interaction_vars, collapse=":"))
    }
  } else if (grepl("^s\\(", term)) {
    plot_title <- paste("CDI Plot for", term)
  } else {
    plot_title <- paste("CDI Plot for", term)
  }
  
  # Create plots
  # 1. Coefficient plot
  coeff_plot <- ggplot(coeffs, aes(x = term_int, y = exp(coeff))) +
    geom_vline(xintercept = unique(coeffs$term_int), color = "grey90") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), width = 0.2, color = obj$colour) +
    geom_point(size = 3, color = obj$colour) +
    scale_x_continuous(breaks = coeffs$term_int, labels = levels(coeffs$term), position = "top") +
    scale_y_log10() +
    labs(x = "", y = "Coefficient") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0),
      panel.grid.minor = element_blank()
    )
  
  # 2. Distribution plot
  dist_plot <- ggplot(distrs, aes(x = term_int, y = focus_int)) +
    geom_vline(xintercept = unique(distrs$term_int), color = "grey90") +
    geom_hline(yintercept = unique(distrs$focus_int), color = "grey90") +
    geom_point(aes(size = prop), color = obj$colour, alpha = 0.7) +
    scale_size_area(max_size = 10) +
    scale_x_continuous(breaks = unique(distrs$term_int), labels = levels(distrs$term)) +
    scale_y_continuous(breaks = unique(distrs$focus_int), labels = levels(distrs$focus)) +
    labs(
      x = variable,
      y = obj$focus
    ) +
    guides(size = "none") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
  
  # 3. Influence plot
  infl_plot <- ggplot(influence_data, aes(x = influence, y = level_int)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = unique(influence_data$level_int), color = "grey90") +
    geom_point(size = 3, color = obj$colour) +
    geom_path(color = obj$colour, group = 1) +
    scale_y_continuous(breaks = influence_data$level_int, labels = influence_data$level, position = "right") +
    labs(x = "Influence", y = "") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank()
    )
  
  # Combine plots
  # Define layout for combined plot - coefficient on top, distribution and influence below
  layout <- "
  AAAAAA
  BBBBCC
  BBBBCC
  "
  
  combined_plot <- coeff_plot + dist_plot + infl_plot + 
    plot_layout(design = layout) +
    plot_annotation(title = plot_title)
  
  return(combined_plot)
}

