#' A package for generating step plots, influence plots, CDI plots, and influence metrics for linear models and GAMs
#'
#' The concept of influence in generalised linear models is described in
#' Bentley, N., Kendrick, T. H., Starr, P. J., & Breen, P. A. (2011). Influence plots and metrics: tools for better understanding fisheries catch-per-unit-effort standardisations.
#' This package provides an implementation of the plots and metrics described in that paper. ICES Journal of Marine Science, doi:10.1093/icesjms/fsr174.
#'
#' This package works for \code{glm} and \code{gam} models with log transformed dependent variables. These
#' are the type of models commonly used for one part of the delta-lognormal approach to catch-per-unit-effort (CPUE)
#' standardisation. e.g.\code{model = gam(log(catch)~s(year)+ti(month,vessel)+effort)} or with interactions
#' \code{model = gam(log(catch)~s(year)+ti(month,year)+effort)}
#'
#' @docType package
#' @name influ-package
#' @aliases influ
#' @title Influence in linear models and GAMs with interaction support
#' @author Nokome Bentley (Original), Alistair Dunn (Modified to support GAMs and ggplot2)

library(ggplot2)
library(dplyr)
library(tidyr)
library(mgcv)
library(gridExtra)
library(patchwork)
library(stringr)

#' Create a new Influence object.
#'
#' A new Influence object needs to be created for each combination of a model and focus term.
#' For example, you might compare the influence of the same term in two separate models:
#' \code{
#'    influ1 = influence(model1, 'year')
#'    influ2 = influence(model2, 'year')
#' }
#'
#' @param model The model for which the influence of terms will be determined
#' @param focus The focus term for which influence is to be calculated
#' @param data The data which was used to fit the model (required for some types of models that do not store the data)
#' @param response The response term in the model
#' @param colour The color to use for plots
#' @return A new Influence object (S3)
#' @export
influence <- function(model, focus = NULL, data = NULL, response = NULL, colour = "steelblue") {
  
  # Create a new influence object
  obj <- structure(list(
    model = model,
    data = data,
    response = response,
    focus = focus,
    colour = colour,
    terms = NULL,
    labels = NULL,
    orders = NULL,
    indices = NULL,
    summary = NULL,
    preds = NULL,
    influences = NULL,
    interactions = NULL  # New field to store interaction information
  ), class = "influence")
  
  # Initialize the object
  obj <- initialize_influence(obj)
  
  return(obj)
}

#' Extract coefficients for a term
#'
#' @param obj The influence object
#' @param model The model to extract coefficients from
#' @param term The term to extract coefficients for
#' @return A vector of coefficients
#' @export
get_coeffs <- function(obj, model = obj$model, term = obj$focus) {
  UseMethod("get_coeffs", model)
}

#' Extract coefficients for a term from a GLM
#'
#' @param obj The influence object
#' @param model The GLM model
#' @param term The term to extract coefficients for
#' @return A vector of coefficients
#' @export
get_coeffs.glm <- function(obj, model = obj$model, term = obj$focus) {
  # Special handling for interaction terms
  if (term %in% names(obj$interactions)) {
    interaction_info <- obj$interactions[[term]]
    if ("year" %in% interaction_info$vars || obj$focus %in% interaction_info$vars) {
      # For interactions involving year or focus variable, use predicted effects
      pred <- get_term_effects(obj, model, term)
      return(pred)
    }
  }
  
  # Standard handling for main effects
  coeffs <- coefficients(model)
  rows <- startsWith(names(coeffs), term)
  c(0, coeffs[rows])
}

#' Default method for get_coeffs
#'
#' @param obj The influence object
#' @param model The model
#' @param term The term to extract coefficients for
#' @return A vector of coefficients
#' @export
get_coeffs.default <- function(obj, model = obj$model, term = obj$focus) {
  coeffs <- coefficients(model)
  rows <- startsWith(names(coeffs), term)
  c(0, coeffs[rows])
}

#' Extract standard errors for coefficients from a GLM
#'
#' @param obj The influence object
#' @param model The GLM model
#' @param term The term to extract SEs for
#' @return Standard errors for the coefficients
#' @export
get_ses.glm <- function(obj, model = obj$model, term = obj$focus) {
  # Special handling for interaction terms
  if (term %in% names(obj$interactions)) {
    interaction_info <- obj$interactions[[term]]
    if ("year" %in% interaction_info$vars || obj$focus %in% interaction_info$vars) {
      # For interactions involving year or focus variable, use predicted standard errors
      pred_data <- create_prediction_grid(obj, term)
      preds <- predict(model, newdata = pred_data, type = "terms", se.fit = TRUE)
      
      # Extract SEs for the term
      if (term %in% colnames(preds$se.fit)) {
        ses <- preds$se.fit[, term]
      } else {
        # For interactions, aggregate SEs
        ses <- rep(mean(diag(vcov(model))), length(get_coeffs(obj, model, term)) - 1)
      }
      
      # Return SEs
      return(ses)
    }
  }
  
  # Standard handling for main effects
  V <- summary(model)$cov.scaled
  rows <- startsWith(row.names(V), term)
  V <- V[rows, rows]
  n <- sum(rows) + 1
  Q <- matrix(-1 / n, nrow = n, ncol = n - 1)
  Q[-1, ] <- Q[-1, ] + diag(rep(1, n - 1))
  V0 <- (Q %*% V) %*% t(Q)
  se <- sqrt(diag(V0))
  se
}

#' Extract standard errors for coefficients from a GAM
#'
#' @param obj The influence object
#' @param model The GAM model
#' @param term The term to extract SEs for
#' @return Standard errors for the coefficients
#' @export
get_ses.gam <- function(obj, model = obj$model, term = obj$focus) {
  # Check if the term is a smooth term or involved in interactions
  if (any(sapply(model$smooth, function(s) s$label == term)) || term %in% names(obj$interactions)) {
    # For smooth terms or interactions, extract standard errors from predict
    pred_data <- create_prediction_grid(obj, term)
    preds <- predict(model, newdata = pred_data, type = "terms", se.fit = TRUE)
    
    # Extract SEs for the term
    if (term %in% colnames(preds$se.fit)) {
      ses <- unique(preds$se.fit[, term])
    } else {
      # For terms not directly in the output (e.g., parts of interactions)
      # Use a conservative estimate based on the average SE across all terms
      all_ses <- as.vector(preds$se.fit)
      ses <- rep(mean(all_ses, na.rm = TRUE), length(get_coeffs(obj, model, term)) - 1)
    }
    
    return(ses)
  } else {
    # For parametric terms not in interactions
    V <- summary(model)$cov.scaled
    rows <- startsWith(row.names(V), term)
    
    if (sum(rows) > 0) {
      V <- V[rows, rows]
      n <- sum(rows) + 1
      Q <- matrix(-1 / n, nrow = n, ncol = n - 1)
      Q[-1, ] <- Q[-1, ] + diag(rep(1, n - 1))
      V0 <- (Q %*% V) %*% t(Q)
      se <- sqrt(diag(V0))
    } else {
      # Fallback for cases where we can't find the term directly
      se <- rep(mean(diag(V)), length(get_coeffs(obj, model, term)) - 1)
    }
    
    return(se)
  }
}

#' Default method for get_ses
#'
#' @param obj The influence object
#' @param model The model
#' @param term The term to extract SEs for
#' @return Standard errors for the coefficients
#' @export
get_ses.default <- function(obj, model = obj$model, term = obj$focus) {
  if (inherits(model, "survreg")) {
    V <- model$var
  } else {
    V <- vcov(model)
  }
  
  rows <- startsWith(row.names(V), term)
  if (sum(rows) > 0) {
    V <- V[rows, rows]
    n <- sum(rows) + 1
    Q <- matrix(-1 / n, nrow = n, ncol = n - 1)
    Q[-1, ] <- Q[-1, ] + diag(rep(1, n - 1))
    V0 <- (Q %*% V) %*% t(Q)
    se <- sqrt(diag(V0))
  } else {
    # Fallback for cases where we can't find the term directly
    se <- rep(mean(diag(V)), length(get_coeffs(obj, model, term)) - 1)
  }
  
  return(se)
}

#' Extract effects for a term
#'
#' @param obj The influence object
#' @param model The model
#' @param term The term to extract effects for
#' @return Effects for the term
#' @export
get_effects <- function(obj, model = obj$model, term = obj$focus) {
  coeffs <- get_coeffs(obj, model, term)
  exp(coeffs - mean(coeffs))
}

#' Print method for influence objects
#'
#' @param x The influence object
#' @param ... Additional arguments
#' @export
print.influence <- function(x, ...) {
  cat("Influence object for model:\n")
  cat("  Response:", x$response, "\n")
  cat("  Focus term:", x$focus, "\n")
  cat("  Terms:", paste(x$terms, collapse=", "), "\n")
  
  # Print interaction information if available
  if (length(x$interactions) > 0) {
    cat("\nInteractions:\n")
    for (term in names(x$interactions)) {
      vars <- paste(x$interactions[[term]]$vars, collapse=", ")
      cat("  ", term, ":", vars, "\n")
    }
  }
  
  if (!is.null(x$summary)) {
    cat("\nSummary statistics:\n")
    print(x$summary)
  } else {
    cat("\nCall calc() to generate summary statistics and enable plotting.\n")
  }
}

#' Standardization plot
#'
#' This plot compares the standardized and unstandardized indices from a model
#'
#' @param obj The influence object
#' @param ... Additional arguments passed to ggplot
#' @return A ggplot object
#' @export
stan_plot <- function(obj, ...) {
  if (is.null(obj$indices)) {
    stop("Call calc() before plotting")
  }
  
  # Prepare data for plotting
  plot_data <- obj$indices %>%
    dplyr::select(level, unstan, stan, stanLower, stanUpper) %>%
    mutate(level_num = as.integer(level))
  
  # Create the plot
  ggplot(plot_data, aes(x = level_num)) +
    geom_errorbar(aes(ymin = stanLower, ymax = stanUpper), width = 0.2) +
    geom_line(aes(y = unstan, group = 1, color = "Unstandardized"), size = 1) +
    geom_point(aes(y = unstan, color = "Unstandardized"), size = 3) +
    geom_line(aes(y = stan, group = 1, color = "Standardized"), size = 1) +
    geom_point(aes(y = stan, color = "Standardized"), size = 3) +
    scale_x_continuous(breaks = plot_data$level_num, labels = plot_data$level) +
    scale_color_manual(values = c("Unstandardized" = "darkgrey", "Standardized" = obj$colour)) +
    labs(
      x = obj$labels[[obj$focus]],
      y = "Index",
      color = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

#' Step plot
#'
#' A plot of the standardized indices as each explanatory variable is added to the model
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

#' Step and influence plots side by side
#'
#' @param obj The influence object
#' @param ... Additional arguments passed to plotting functions
#' @return A combined plot
#' @export
step_and_influ_plot <- function(obj, ...) {
  # Create both plots
  step_p <- step_plot(obj, panels = TRUE, ...)
  influ_p <- influ_plot(obj, panels = TRUE, ...)
  
  # Combine them side by side
  combined_plot <- step_p + influ_p + plot_layout(ncol = 2)
  return(combined_plot)
}

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

#' Plot a CDI plot for each of the terms in the model
#'
#' @param obj The influence object
#' @param ... Additional arguments passed to cdi_plot
#' @return A list of plots
#' @export
cdi_plot_all <- function(obj, ...) {
  plots <- list()
  
  for (term in obj$terms) {
    if (term != obj$focus) {
      plots[[term]] <- cdi_plot(obj, term, ...)
    }
  }
  
  return(plots)
}

# The initialize_influence function needs to handle 'year' as a factor
initialize_influence <- function(obj) {
  # If no data was supplied
  if (is.null(obj$data)) {
    # See if the model has a data variable
    obj$data <- obj$model$data
    if (is.null(obj$data)) {
      # See if it is available as a model attribute
      obj$data <- attr(obj$model, "data")
    }
    if (is.null(obj$data)) {
      # For GAMs, use model$model
      if (inherits(obj$model, "gam")) {
        obj$data <- obj$model$model
      } else {
        stop("You need to provide the data that was fitted to for this type of model e.g. influence(model, data=mydata, ...)")
      }
    }
  }
  
  # Get term labels based on model type
  if (inherits(obj$model, "gam")) {
    # For GAMs, handle smooth terms differently
    parameterization <- obj$model$pterms
    obj$terms <- attr(parameterization, "term.labels")
    
    # Add smooth terms from the model
    smooth_terms <- names(obj$model$smooth)
    if (length(smooth_terms) > 0) {
      # Store interaction information
      obj$interactions <- list()
      
      for (i in 1:length(obj$model$smooth)) {
        smooth_label <- obj$model$smooth[[i]]$label
        obj$terms <- c(obj$terms, smooth_label)
        
        # Check for interactions with year or focus
        if (inherits(obj$model$smooth[[i]], "tensor.smooth") || 
            inherits(obj$model$smooth[[i]], "tp") || 
            inherits(obj$model$smooth[[i]], "ti")) {
          
          # Get variables involved in this smooth
          vars <- obj$model$smooth[[i]]$term
          
          # Check if this is an interaction involving the focus term or year
          if (length(vars) > 1) {
            # Store interaction information
            obj$interactions[[smooth_label]] <- list(
              type = class(obj$model$smooth[[i]])[1],
              vars = vars,
              by = obj$model$smooth[[i]]$by,
              label = smooth_label
            )
          }
        }
      }
    }
  } else {
    # For GLMs and other models
    obj$terms <- attr(obj$model$terms, "term.labels")
    
    # Check for interactions with year or focus in GLM terms
    obj$interactions <- list()
    for (term in obj$terms) {
      if (grepl(":", term, fixed = TRUE)) {
        # This is an interaction term
        parts <- strsplit(term, ":", fixed = TRUE)[[1]]
        obj$interactions[[term]] <- list(
          type = "parametric",
          vars = parts,
          label = term
        )
      }
    }
  }
  
  # Set response and focus if not provided
  if (is.null(obj$response)) {
    obj$response <- as.character(formula(obj$model)[[2]])
  }
  
  if (is.null(obj$focus)) {
    obj$focus <- obj$terms[1]
  }
  
  # Initialize labels and orders
  obj$labels <- setNames(as.list(obj$terms), obj$terms)
  obj$orders <- setNames(as.list(rep("asis", length(obj$terms))), obj$terms)
  
  return(obj)
}

# Update create_prediction_grid to properly handle factor year
create_prediction_grid <- function(obj, term) {
  # Initialize an empty list to store the grid variables
  grid_vars <- list()
  
  # Always include the focus variable - ensuring factor handling
  focus_var <- obj$data[[obj$focus]]
  if (is.factor(focus_var)) {
    grid_vars[[obj$focus]] <- levels(focus_var)
  } else if (is.numeric(focus_var)) {
    grid_vars[[obj$focus]] <- seq(min(focus_var), max(focus_var), length.out = 20)
  } else {
    grid_vars[[obj$focus]] <- unique(focus_var)
  }
  
  # If the term is an interaction, include all variables involved
  if (term %in% names(obj$interactions)) {
    interaction_vars <- obj$interactions[[term]]$vars
    for (var in interaction_vars) {
      if (var != obj$focus && !(var %in% names(grid_vars))) {
        var_data <- obj$data[[var]]
        if (is.factor(var_data)) {
          grid_vars[[var]] <- levels(var_data)
        } else if (is.numeric(var_data)) {
          grid_vars[[var]] <- seq(min(var_data), max(var_data), length.out = 10)
        } else {
          grid_vars[[var]] <- unique(var_data)
        }
      }
    }
  }
  
  # For all other terms involved in interactions with the focus or year
  if (term == obj$focus || term == "year") {
    for (interaction_name in names(obj$interactions)) {
      interaction_info <- obj$interactions[[interaction_name]]
      if (obj$focus %in% interaction_info$vars || "year" %in% interaction_info$vars) {
        for (var in interaction_info$vars) {
          if (var != obj$focus && var != "year" && !(var %in% names(grid_vars))) {
            var_data <- obj$data[[var]]
            if (is.factor(var_data)) {
              grid_vars[[var]] <- levels(var_data)
            } else if (is.numeric(var_data)) {
              grid_vars[[var]] <- seq(min(var_data), max(var_data), length.out = 10)
            } else {
              grid_vars[[var]] <- unique(var_data)
            }
          }
        }
      }
    }
  }
  
  # Handle year variable separately if it's not the focus
  if ("year" %in% names(obj$data) && !"year" %in% names(grid_vars) && "year" != obj$focus) {
    year_data <- obj$data[["year"]]
    if (is.factor(year_data)) {
      grid_vars[["year"]] <- levels(year_data)
    } else if (is.numeric(year_data)) {
      grid_vars[["year"]] <- seq(min(year_data), max(year_data), length.out = 10)
    } else {
      grid_vars[["year"]] <- unique(year_data)
    }
  }
  
  # Add all other variables needed for prediction
  all_vars <- names(obj$data)
  for (var in all_vars) {
    if (!(var %in% names(grid_vars)) && var != obj$response) {
      var_data <- obj$data[[var]]
      # Use median/mode values for variables not directly involved
      if (is.factor(var_data)) {
        # Use the most common level (mode)
        grid_vars[[var]] <- names(sort(table(var_data), decreasing = TRUE)[1])
      } else if (is.numeric(var_data)) {
        # Use median for numeric variables
        grid_vars[[var]] <- median(var_data)
      } else {
        # Use the first value for other types
        grid_vars[[var]] <- var_data[1]
      }
    }
  }
  
  # Create the grid using expand.grid
  pred_grid <- do.call(expand.grid, grid_vars)
  
  # Ensure all factors are preserved as factors with proper levels
  for (var in names(pred_grid)) {
    if (is.factor(obj$data[[var]])) {
      pred_grid[[var]] <- factor(pred_grid[[var]], levels = levels(obj$data[[var]]))
    }
  }
  
  return(pred_grid)
}

# Update get_coeffs for GAMs to handle year as factor in smooths
get_coeffs.gam <- function(obj, model = obj$model, term = obj$focus) {
  # Check if the term is a smooth term
  if (any(sapply(model$smooth, function(s) s$label == term))) {
    # For smooth terms, extract coefficients from predict
    pred <- get_term_effects(obj, model, term)
    return(pred)
  } else {
    # For parametric terms
    coeffs <- coefficients(model)
    rows <- startsWith(names(coeffs), term)
    if (sum(rows) > 0) {
      c(0, coeffs[rows])
    } else {
      # Handle case where term might be part of an interaction but not a main effect
      pred <- get_term_effects(obj, model, term)
      return(pred)
    }
  }
}

# Updated get_term_effects to properly handle year as a factor
get_term_effects <- function(obj, model, term) {
  # Check if we're dealing with a focus term that's involved in interactions
  is_interaction <- term %in% names(obj$interactions)
  has_interactions_with_focus <- any(sapply(obj$interactions, function(x) obj$focus %in% x$vars))
  has_interactions_with_year <- any(sapply(obj$interactions, function(x) "year" %in% x$vars))
  
  # Create a prediction grid
  if (is_interaction || (term == obj$focus && has_interactions_with_focus) || 
      (term == "year" && has_interactions_with_year)) {
    # For interactions or focus terms with interactions
    pred_data <- create_prediction_grid(obj, term)
    
    # Get predictions
    preds <- predict(model, newdata = pred_data, type = "terms", se.fit = TRUE)
    
    # Extract predictions for the term
    if (term %in% colnames(preds$fit)) {
      # Direct term effect
      effects <- preds$fit[, term]
    } else if (is_interaction) {
      # Interaction effect
      interaction_vars <- obj$interactions[[term]]$vars
      
      # For tensor smooths or ti() terms
      if (obj$interactions[[term]]$type %in% c("tensor.smooth", "tp", "ti")) {
        effects <- preds$fit[, term]
      } else {
        # For parametric interactions, we need to aggregate over the non-focus variables
        focus_var <- if(obj$focus %in% interaction_vars) obj$focus else 
                     if("year" %in% interaction_vars) "year" else interaction_vars[1]
        
        # Group by focus variable and average other effects
        agg_effects <- aggregate(
          preds$fit[, term], 
          by = list(focus_level = pred_data[[focus_var]]), 
          FUN = mean
        )
        effects <- agg_effects$x
      }
    } else {
      # No direct effect, might be marginal effect from interactions
      # Sum all effects related to the term
      related_terms <- colnames(preds$fit)[grepl(term, colnames(preds$fit))]
      if (length(related_terms) > 0) {
        effects_mat <- preds$fit[, related_terms, drop = FALSE]
        effects <- rowSums(effects_mat)
      } else {
        # Fallback
        effects <- rep(0, nrow(pred_data))
      }
    }
    
    # Aggregate effects by focus variable
    if (term == obj$focus || (term == "year" && obj$focus != "year")) {
      # For the focus term or year, return the effects directly
      unique_effects <- unique(effects)
      return(c(0, unique_effects))
    } else {
      # For other terms, return the aggregated effects
      unique_effects <- unique(effects)
      return(c(0, unique_effects))
    }
  } else {
    # For non-interaction terms or terms not involving the focus or year
    if (inherits(model, "gam") && any(sapply(model$smooth, function(s) s$label == term))) {
      # For smooth terms
      newdata <- obj$data
      term_effect <- predict(model, newdata, type = "terms", terms = term)
      # Return unique levels with 0 prepended
      unique_levels <- unique(term_effect)
      return(c(0, unique_levels))
    } else {
      # For parametric terms
      coeffs <- coefficients(model)
      rows <- startsWith(names(coeffs), term)
      return(c(0, coeffs[rows]))
    }
  }
}

# Update get_ses.gam for proper handling of year as factor
get_ses.gam <- function(obj, model = obj$model, term = obj$focus) {
  # Check if the term is a smooth term or involved in interactions
  if (any(sapply(model$smooth, function(s) s$label == term)) || term %in% names(obj$interactions) || term == "year") {
    # For smooth terms or interactions, extract standard errors from predict
    pred_data <- create_prediction_grid(obj, term)
    preds <- predict(model, newdata = pred_data, type = "terms", se.fit = TRUE)
    
    # Extract SEs for the term
    if (term %in% colnames(preds$se.fit)) {
      ses <- unique(preds$se.fit[, term])
    } else {
      # For terms not directly in the output (e.g., parts of interactions)
      # Use a conservative estimate based on the average SE across all terms
      all_ses <- as.vector(preds$se.fit)
      ses <- rep(mean(all_ses, na.rm = TRUE), length(get_coeffs(obj, model, term)) - 1)
    }
    
    return(ses)
  } else {
    # For parametric terms not in interactions
    V <- summary(model)$cov.scaled
    rows <- startsWith(row.names(V), term)
    
    if (sum(rows) > 0) {
      V <- V[rows, rows]
      n <- sum(rows) + 1
      Q <- matrix(-1 / n, nrow = n, ncol = n - 1)
      Q[-1, ] <- Q[-1, ] + diag(rep(1, n - 1))
      V0 <- (Q %*% V) %*% t(Q)
      se <- sqrt(diag(V0))
    } else {
      # Fallback for cases where we can't find the term directly
      se <- rep(mean(diag(V)), length(get_coeffs(obj, model, term)) - 1)
    }
    
    return(se)
  }
}

# Update the calc function to better handle year as a factor
calc <- function(obj, islog = NULL) {
  # Get observed values
  if (inherits(obj$model, "gam")) {
    observed <- obj$model$model[[obj$response]]
  } else {
    observed <- obj$model$model[, obj$response]
  }
  
  if (inherits(observed, "Surv")) {
    observed <- as.numeric(observed[, 1])
  }
  
  if (is.null(islog)) {
    logged <- startsWith(obj$response, "log(")
  } else {
    logged <- islog
  }
  
  # Extract focus variable levels
  if (inherits(obj$model, "gam")) {
    focus_var <- obj$model$model[[obj$focus]]
  } else {
    focus_var <- obj$model$model[, obj$focus]
  }
  
  # Ensure focus_var is a factor
  if (!is.factor(focus_var)) {
    focus_var <- factor(focus_var)
  }
  
  # Create indices dataframe
  obj$indices <- data.frame(level = levels(focus_var))
  
  # Add unstandardized index
  if (logged || sum(observed <= 0) == 0) {
    if (logged) {
      log_observed <- observed
    } else {
      log_observed <- log(observed)
    }
    # Calculate geometric mean
    unstan_df <- aggregate(list(unstan = log_observed), list(level = focus_var), mean)
    obj$indices <- merge(obj$indices, unstan_df)
    # Turn into relative index
    obj$indices$unstan <- with(obj$indices, exp(unstan - mean(unstan)))
  } else {
    # Use arithmetic mean instead
    unstan_df <- aggregate(list(unstan = observed), list(level = focus_var), mean)
    obj$indices <- merge(obj$indices, unstan_df)
    # Turn into relative index
    obj$indices$unstan <- with(obj$indices, unstan / mean(unstan))
  }
  
  # Add standardized index
  coeffs <- get_coeffs(obj)
  ses <- get_ses(obj)
  base <- mean(coeffs)
  obj$indices$stan <- exp(coeffs - base)
  obj$indices$stanLower <- exp(coeffs - base - 2 * ses)
  obj$indices$stanUpper <- exp(coeffs - base + 2 * ses)
  
  # Create models with terms successively added
  obj$summary <- NULL
  
  # Calculate influence statistics
  for (termCount in 0:length(obj$terms)) {
    if (termCount > 0) {
      term <- obj$terms[termCount]
      
      # Handle models with interactions differently
      if (inherits(obj$model, "gam") && any(sapply(obj$model$smooth, function(s) s$label == term))) {
        # For smooth terms in GAMs
        terms_to_include <- obj$terms[1:termCount]
        smooth_terms <- sapply(obj$model$smooth, function(s) s$label)
        smooth_to_include <- smooth_terms[smooth_terms %in% terms_to_include]
        
        # Build formula for GAM, preserving interactions
        formula_terms <- character(0)
        for (t in terms_to_include) {
          if (t %in% smooth_to_include) {
            # Find the original specification of the smooth term
            for (s in obj$model$smooth) {
              if (s$label == t) {
                # Handle factor year in smooth terms
                if ("year" %in% s$term && is.factor(obj$data[["year"]])) {
                  # For factor year in smooths, convert to parametric
                  if (inherits(s, "tensor.smooth") || inherits(s, "ti")) {
                    # For interactions with year as factor, modify the term
                    vars <- s$term
                    year_idx <- which(vars == "year")
                    other_vars <- vars[-year_idx]
                    
                    if (length(other_vars) > 0) {
                      # For interactions with other variables
                      other_var_terms <- paste0("s(", other_vars, ")")
                      year_term <- "year"
                      # Create interaction between year and smooth
                      formula_terms <- c(formula_terms, 
                                         other_var_terms,
                                         year_term,
                                         paste0(other_var_terms, ":year"))
                    } else {
                      # Just year by itself
                      formula_terms <- c(formula_terms, "year")
                    }
                  } else {
                    # Single term with year as factor
                    formula_terms <- c(formula_terms, "year")
                  }
                } else {
                  # Regular smooth terms without factor year
                  if (inherits(s, "tensor.smooth") || inherits(s, "ti")) {
                    # For tensor interaction smooth terms (ti)
                    vars <- paste(s$term, collapse = ",")
                    formula_terms <- c(formula_terms, paste0("ti(", vars, ")"))
                  } else {
                    # For regular smooth terms (s)
                    formula_terms <- c(formula_terms, paste0("s(", s$term, ")"))
                  }
                }
                break
              }
            }
          } else if (t %in% names(obj$interactions) && obj$interactions[[t]]$type == "parametric") {
            # For parametric interactions
            formula_terms <- c(formula_terms, t)
          } else {
            # For regular parametric terms
            formula_terms <- c(formula_terms, t)
          }
        }
        
        # Create the updated formula
        formula_rhs <- paste(formula_terms, collapse = "+")
        formula_str <- paste(obj$response, "~", formula_rhs)
        new_formula <- as.formula(formula_str)
        
        # Update model with new formula
        model <- update(obj$model, formula = new_formula)
      } else {
        # For GLMs and parametric terms, preserving interactions
        terms_to_include <- obj$terms[1:termCount]
        
        # Check for interaction terms involving any of the included terms
        interaction_terms <- names(obj$interactions)[sapply(obj$interactions, function(i) {
          any(i$vars %in% terms_to_include)
        })]
        
        # Only include interaction terms if all their variables are in terms_to_include
        interaction_terms_to_include <- character(0)
        for (i_term in interaction_terms) {
          if (all(obj$interactions[[i_term]]$vars %in% terms_to_include)) {
            interaction_terms_to_include <- c(interaction_terms_to_include, i_term)
          }
        }
        
        # Combine main effects and allowed interaction terms
        all_terms_to_include <- unique(c(
          terms_to_include[!terms_to_include %in% names(obj$interactions)],
          interaction_terms_to_include
        ))
        
        # Create formula
        formula_rhs <- paste(all_terms_to_include, collapse = "+")
        formula_str <- paste(obj$response, "~", formula_rhs)
        new_formula <- as.formula(formula_str)
        
        # Update model
        model <- update(obj$model, formula = new_formula)
      }
      
      # Get index for this model
      index <- get_effects(obj, model)
      
      # Add column to indices
      obj$indices <- cbind(obj$indices, index)
      
      # Give index column the right hand side of formula as name
      names(obj$indices)[ncol(obj$indices)] <- if (termCount == 1) term else paste("+", term)
    } else {
      term <- "intercept"
      
      # Create intercept-only model
      if (inherits(obj$model, "gam")) {
        formula_str <- paste(obj$response, "~ 1")
        model <- gam(as.formula(formula_str), data = obj$data, family = obj$model$family)
      } else {
        model <- update(obj$model, . ~ 1)
      }
    }
    
    # Extract model statistics
    if (inherits(model, "survreg")) {
      logLike <- model$loglik[2]
      fitted <- predict(model, type = "response")
    } else {
      logLike <- as.numeric(logLik(model))
      fitted <- fitted(model)
    }
    
    # Calculate R-squared values
    if (termCount == 0) {
      r2 <- 0
    } else {
      if (!logged) {
        r2 <- cor(log(observed), log(fitted))^2
      } else {
        r2 <- cor(observed, fitted)^2
      }
    }
    
    # Deviance pseudo-R2
    if (inherits(model, c("glm", "gam"))) {
      r2Dev <- (model$null.deviance - model$deviance) / model$null.deviance
    } else {
      r2Dev <- NA
    }
    
    # Negelkerke pseudo-R2
    if (termCount == 0) {
      logLikeInterceptOnly <- logLike
    }
    n <- length(observed)
    r2Negel <- (1 - exp((logLikeInterceptOnly - logLike) * (2 / n))) / (1 - exp(logLikeInterceptOnly * (2 / n)))
    
    # Add to summary dataframe
    obj$summary <- rbind(obj$summary, data.frame(
      term = term,
      k = length(coef(model)),
      logLike = logLike,
      aic = extractAIC(model)[2],
      r2 = r2,
      r2Dev = r2Dev,
      r2Negel = r2Negel
    ))
  }
  
  # R2 values presented as differences
  obj$summary <- within(obj$summary, {
    k <- c(1, diff(k))
    r2 <- c(NA, diff(r2))
    r2Dev <- c(NA, diff(r2Dev))
    r2Negel <- c(NA, diff(r2Negel))
  })
  
  # Calculate predicted values and handle interactions
  if (inherits(obj$model, "gam")) {
    # For GAMs, get predictions for each term separately
    preds_fit <- list()
    preds_se <- list()
    
    # Handle parametric terms first
    if (length(obj$model$pterms) > 0) {
      param_terms <- attr(obj$model$pterms, "term.labels")
      for (term in param_terms) {
        if (term != "") {  # Skip empty terms
          # Check if this is part of an interaction
          if (term %in% names(obj$interactions)) {
            # For interaction terms, we need a prediction grid
            pred_data <- create_prediction_grid(obj, term)
            pred <- predict(obj$model, newdata = pred_data, type = "terms", terms = term, se.fit = TRUE)
            
            # Map back to the original data points
            # This is a simplification - in practice we'd need to interpolate or match values
            all_preds <- predict(obj$model, type = "terms", terms = term, se.fit = TRUE)
            preds_fit[[term]] <- all_preds$fit
            preds_se[[term]] <- all_preds$se.fit
          } else {
            # For regular terms
            pred <- predict(obj$model, type = "terms", terms = term, se.fit = TRUE)
            preds_fit[[term]] <- pred$fit
            preds_se[[term]] <- pred$se.fit
          }
        }
      }
    }
    
    # Handle smooth terms
    for (i in seq_along(obj$model$smooth)) {
      term <- obj$model$smooth[[i]]$label
      
      # Check if this is an interaction term involving the focus or year
      if (term %in% names(obj$interactions) && 
          (obj$focus %in% obj$interactions[[term]]$vars || "year" %in% obj$interactions[[term]]$vars)) {
        # For interaction terms involving focus or year
        pred_data <- create_prediction_grid(obj, term)
        pred <- predict(obj$model, newdata = pred_data, type = "terms", terms = term, se.fit = TRUE)
        
        # Map back to original data - simplified approach
        all_preds <- predict(obj$model, type = "terms", terms = term, se.fit = TRUE)
        preds_fit[[term]] <- all_preds$fit
        preds_se[[term]] <- all_preds$se.fit
      } else {
        # For regular smooth terms
        pred <- predict(obj$model, type = "terms", terms = term, se.fit = TRUE)
        preds_fit[[term]] <- pred$fit
        preds_se[[term]] <- pred$se.fit
      }
    }
    
    # Create dataframe with all predictions
    preds_df <- data.frame(obj$data)
    for (term in names(preds_fit)) {
      preds_df[[paste0("fit.", term)]] <- preds_fit[[term]]
      preds_df[[paste0("se.fit.", term)]] <- preds_se[[term]]
    }
    obj$preds <- preds_df
  } else {
    # For GLMs, handle interactions specifically
    has_interactions <- length(obj$interactions) > 0
    
    if (has_interactions) {
      # Get standard predictions for all terms
      preds <- predict(obj$model, type = "terms", se.fit = TRUE)
      fit <- as.data.frame(preds$fit)
      se.fit <- as.data.frame(preds$se.fit)
      
      # Process interactions involving the focus or year
      for (term in names(obj$interactions)) {
        interaction_info <- obj$interactions[[term]]
        if (obj$focus %in% interaction_info$vars || "year" %in% interaction_info$vars) {
          # For interactions involving focus or year, we need special handling
          if (term %in% colnames(fit)) {
            # The term exists in the predictions, no special handling needed
            next
          }
          
          # Create grid for this interaction
          pred_data <- create_prediction_grid(obj, term)
          
          # Get predictions on the grid
          pred_grid <- predict(obj$model, newdata = pred_data, type = "terms", se.fit = TRUE)
          
          # Map back to original data points - simplified approach
          # In a real implementation, we'd need to match or interpolate values
          if (term %in% colnames(pred_grid$fit)) {
            # Add as new column in the fit and se.fit dataframes
            grid_values <- pred_grid$fit[, term]
            grid_se <- pred_grid$se.fit[, term]
            
            # Map values to original data by matching focus or year levels
            mapped_values <- numeric(nrow(obj$data))
            mapped_se <- numeric(nrow(obj$data))
            
            # Determine which variable to use for mapping
            map_var <- if (obj$focus %in% interaction_info$vars) obj$focus else "year"
            
            for (level in levels(factor(obj$data[[map_var]]))) {
              level_rows <- obj$data[[map_var]] == level
              grid_rows <- pred_data[[map_var]] == level
              
              if (any(level_rows) && any(grid_rows)) {
                # Take the mean value for this level
                mapped_values[level_rows] <- mean(grid_values[grid_rows])
                mapped_se[level_rows] <- mean(grid_se[grid_rows])
              }
            }
            
            # Add to fit and se.fit
            fit[[term]] <- mapped_values
            se.fit[[term]] <- mapped_se
          }
        }
      }
      
      # Combine into preds dataframe
      preds <- cbind(fit, se.fit)
      names(preds) <- c(paste0("fit.", names(fit)), paste0("se.fit.", names(se.fit)))
      obj$preds <- cbind(obj$data, preds)
    } else {
      # Standard GLM without interactions
      preds <- predict(obj$model, type = "terms", se.fit = TRUE)
      fit <- as.data.frame(preds$fit)
      se.fit <- as.data.frame(preds$se.fit)
      preds <- cbind(fit, se.fit)
      names(preds) <- c(paste0("fit.", names(fit)), paste0("se.fit.", names(se.fit)))
      obj$preds <- cbind(obj$data, preds)
    }
  }
  
  # Calculate influences and statistics for each term
  focus_var <- obj$data[[obj$focus]]
  if (!is.factor(focus_var)) {
    focus_var <- factor(focus_var)
  }
  
  obj$influences <- data.frame(level = levels(focus_var))
  overall <- c(NA, NA)  # NAs for null model and for focus term
  trend <- c(NA, NA)
  
  # Handle influences for interactions with year or focus
  for (term in obj$terms) {
    if (term != obj$focus) {
      # Check if this is an interaction term involving focus or year
      if (term %in% names(obj$interactions) && 
          (obj$focus %in% obj$interactions[[term]]$vars || "year" %in% obj$interactions[[term]]$vars)) {
        
        # For interaction terms, we need to use the prediction grid
        pred_data <- create_prediction_grid(obj, term)
        pred <- predict(obj$model, newdata = pred_data, type = "terms", terms = term)
        
        # Determine which variable to use for aggregation (focus or year)
        agg_var <- if (obj$focus %in% obj$interactions[[term]]$vars) obj$focus else "year"
        
        # Aggregate by focus or year level
        infl <- aggregate(
          list(value = pred),
          by = list(level = pred_data[[agg_var]]),
          mean
        )
        
        # If we aggregated by year but focus is not year, map to focus levels
        if (agg_var == "year" && obj$focus != "year") {
          # Create a mapping from year to focus levels based on the original data
          year_focus_map <- aggregate(
            list(focus_level = obj$data[[obj$focus]]),
            by = list(year_level = obj$data[["year"]]),
            function(x) names(sort(table(x), decreasing = TRUE))[1]  # Most common focus level for each year
          )
          
          # Apply the mapping to get focus levels
          infl$year_level <- infl$level
          infl$level <- sapply(infl$year_level, function(y) {
            idx <- which(year_focus_map$year_level == y)
            if (length(idx) > 0) year_focus_map$focus_level[idx[1]] else NA
          })
          
          # Aggregate again by focus level
          infl <- aggregate(
            list(value = infl$value),
            by = list(level = infl$level),
            mean
          )
        }
        
        overall <- c(overall, with(infl, exp(mean(abs(value))) - 1))
        trend <- c(trend, with(infl, exp(cov(1:length(value), value) / var(1:length(value))) - 1))
        names(infl) <- c("level", term)
        obj$influences <- merge(obj$influences, infl, all.x = TRUE, by = "level")
      } else {
        # For regular terms, use standard approach
        term_col <- paste0("fit.", term)
        
        if (term_col %in% names(obj$preds)) {
          infl <- aggregate(
            list(value = obj$preds[[term_col]]),
            list(level = obj$preds[[obj$focus]]),
            mean
          )
          overall <- c(overall, with(infl, exp(mean(abs(value))) - 1))
          trend <- c(trend, with(infl, exp(cov(1:length(value), value) / var(1:length(value))) - 1))
          names(infl) <- c("level", term)
          obj$influences <- merge(obj$influences, infl, all.x = TRUE, by = "level")
        }
      }
    }
  }
  
  
  # Insert statistics into summary table
  obj$summary$overall <- overall
  obj$summary$trend <- trend
  
  return(obj)
}

# Example usage:
# 
# # Load required packages
# library(mgcv)
# library(ggplot2)
# library(dplyr)
# 
# # Create example data with year interactions
# set.seed(123)
# n <- 500
# data <- data.frame(
#   year = factor(rep(2010:2019, each = n/10)),
#   month = factor(rep(1:12, n/12)),
#   vessel = factor(rep(1:5, each = n/5)),
#   area = factor(rep(1:3, each = n/3)),
#   effort = runif(n, 0.5, 10),
#   catch = NA
# )
# 
# # Generate catch data with year-month interaction effects
# year_effect <- c(1.0, 1.2, 1.4, 1.3, 0.9, 0.8, 0.7, 0.8, 0.9, 1.0)
# month_effect <- c(0.7, 0.8, 1.0, 1.2, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6)
# vessel_effect <- c(0.8, 1.0, 1.2, 0.9, 1.1)
# area_effect <- c(0.9, 1.0, 1.1)
# 
# # Create year-month interaction effect matrix
# year_month_interaction <- matrix(1, nrow = length(year_effect), ncol = length(month_effect))
# # Add some interaction patterns (e.g., certain months have different patterns in different years)
# for (y in 1:length(year_effect)) {
#   for (m in 1:length(month_effect)) {
#     # Create a seasonal pattern that shifts over the years
#     year_month_interaction[y, m] <- 1 + 0.3 * sin((y + m/12) * pi/5)
#   }
# }
# 
# for (i in 1:n) {
#   year_idx <- as.integer(data$year[i]) - 2009
#   month_idx <- as.integer(data$month[i])
#   vessel_idx <- as.integer(data$vessel[i])
#   area_idx <- as.integer(data$area[i])
#   
#   # Base catch with effort effect
#   base_catch <- data$effort[i] * 2
#   
#   # Apply effects and interactions
#   interaction_effect <- year_month_interaction[year_idx, month_idx]
#   true_catch <- base_catch * year_effect[year_idx] * month_effect[month_idx] * 
#                vessel_effect[vessel_idx] * area_effect[area_idx] * interaction_effect
#   
#   # Add noise
#   data$catch[i] <- true_catch * exp(rnorm(1, 0, 0.2))
# }
# 
# # Log-transform catch for modeling
# data$log_catch <- log(data$catch)
# 
# # Model with parametric year:month interaction
# glm_model <- glm(log_catch ~ year + month + year:month + vessel + area + log(effort), data = data)
# 
# # Create influence object and calculate statistics
# infl_glm <- influence(glm_model, focus = "year")
# infl_glm <- calc(infl_glm)
# 
# # Create plots
# stan_plot(infl_glm)
# step_plot(infl_glm)
# influ_plot(infl_glm)
# cdi_plot(infl_glm, "year:month")
# 
# # GAM model with tensor product interaction between year and month
# gam_model <- gam(log_catch ~ s(year, k=10) + s(month, k=12) + ti(year, month, k=c(5,5)) + 
#                 s(vessel) + s(area) + log(effort), 
#                 data = data, method = "REML")
# 
# # Create influence object for GAM and calculate statistics
# infl_gam <- influence(gam_model, focus = "year")
# infl_gam <- calc(infl_gam)
# 
# # Create plots for GAM with tensor interactions
# stan_plot(infl_gam)
# step_plot(infl_gam) 
# influ_plot(infl_gam)
# cdi_plot(infl_gam, "ti(year,month)")

# For GLMs with parametric interactions
infl_glm <- influence(gam1, focus = "year")
create_prediction_grid(infl_glm, "year")
infl_glm <- calc(infl_glm)

infl_glm <- calc(infl_glm)
step_plot(infl_glm)

# For GAMs with tensor smooth interactions
infl_gam <- influence(gam3, focus = "year")
infl_gam <- calc(infl_gam)
step_plot(infl_gam)