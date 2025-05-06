# Load libraries (Imports preferred in packages)
library(ggplot2)
library(tidyr)
library(patchwork)
library(RColorBrewer)

# Source helper functions (if they are in separate files)
# source("R/plot-helpers.R") # Assuming helpers are grouped

# -----------------------------------------------------------------------------
# Plotting Method (Rewritten for ggplot2)
# -----------------------------------------------------------------------------

# Define generic if it doesn't exist (it does in base R)
# plot <- function(x, ...) UseMethod("plot")

#' Plot influence_gam Object using ggplot2
#'
#' Generates various diagnostic plots for influence analysis results from a
#' `influence_gam` object using the `ggplot2` package.
#'
#' @param x An object of class `influence_gam` with calculated results.
#' @param y Ignored. Included for compatibility with the generic plot function.
#' @param type The type of plot to generate. Options are:
#'   `"stan"`: Standardization plot (unstandardised vs standardised index).
#'   `"step"`: Step plot showing index changes as terms are added.
#'   `"influ"`: Influence plot showing the effect of each non-focus term across focus levels.
#'   `"cdi"`: Coefficient-Distribution-Influence plot (EXPERIMENTAL for GAMs).
#'   `"step_influ"`: Side-by-side step and influence plots (using patchwork).
#'   `"interaction"`: Plot for visualizing two-way factor interaction terms.
#' @param term For `type = "cdi"` or `type = "interaction"`, the specific model term (character string) to plot.
#' @param panels For `type = "step"` or `type = "influ"`, should plots be panelled (TRUE, using facets) or overlaid (FALSE)?
#' @param labels A named list to override default axis labels (e.g., `list(focus = "Year", term = "Effort")`). Names should match focus variable or term names.
#' @param base_size Base font size for ggplot theme.
#' @param ... Additional arguments (passed to helper functions or underlying methods, currently less used).
#'
#' @return A ggplot object (for single plots) or a patchwork object (for combined plots like cdi, step_influ).
#' @method plot influence_gam
#' @export
plot.influence_gam <- function(x, y = NULL, # Add y for generic compatibility
                               type = c("stan", "step", "influ", "cdi", "step_influ", "interaction"),
                               term = NULL,
                               panels = TRUE,
                               labels = list(),
                               base_size = 11,
                               ...) {
  if (!inherits(x, "influence_gam")) {
    stop("Input 'x' must be of class 'influence_gam'.")
  }
  if (!x$calculated) {
    stop("Calculations not performed. Run 'calculate_influence()' before plotting.")
  }

  # Check if necessary plotting packages are loaded/available
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' needed for plotting.")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Package 'RColorBrewer' needed for plotting.")
  if (type %in% c("step", "influ") && !requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' needed for 'step' and 'influ' plots.")
  if (type %in% c("cdi", "step_influ") && !requireNamespace("patchwork", quietly = TRUE)) stop("Package 'patchwork' needed for 'cdi' and 'step_influ' plots.")


  type <- match.arg(type)
  obj <- x # Use obj internally for clarity
  focus_var <- obj$focus_var
  focus_label <- labels[[focus_var]] %||% focus_var # Use custom label or variable name

  # --- Plotting Logic ---

  # Standardization Plot
  if (type == "stan") {
    if (!all(c("level", "unstan", "stan", "stanLower", "stanUpper") %in% names(obj$indices))) {
      stop("Required columns missing in obj$indices for stan plot.")
    }
    indices <- obj$indices
    indices$level_num <- as.numeric(factor(indices$level)) # Ensure level_num corresponds to factor levels

    p <- ggplot(indices, aes(x = level_num)) +
      geom_point(aes(y = unstan, shape = "Unstandardised", color = "Unstandardised")) +
      geom_line(aes(y = unstan, linetype = "Unstandardised", color = "Unstandardised")) +
      geom_point(aes(y = stan, shape = "Standardised", color = "Standardised")) +
      geom_line(aes(y = stan, linetype = "Standardised", color = "Standardised")) +
      geom_errorbar(aes(ymin = stanLower, ymax = stanUpper, color = "Standardised"), width = 0.2, na.rm = TRUE) + # Add na.rm
      scale_x_continuous(breaks = indices$level_num, labels = indices$level) +
      scale_shape_manual("Index Type", values = c("Unstandardised" = 1, "Standardised" = 16)) +
      scale_linetype_manual("Index Type", values = c("Unstandardised" = "dashed", "Standardised" = "solid")) + # Diff linetypes
      scale_color_brewer("Index Type", palette = "Set1") + # Use color palette
      labs(title = "Standardization Plot", x = focus_label, y = "Index (Relative Scale)") +
      theme(legend.position = "bottom")
    return(p)
  }

  # Step Plot
  else if (type == "step") {
    step_data_long <- prepare_step_data(obj$indices, focus_var)
    if (is.null(step_data_long)) stop("No step index columns found in obj$indices.")

    p <- ggplot(step_data_long, aes(x = level_num, y = index_value)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
      scale_x_continuous(breaks = unique(step_data_long$level_num), labels = levels(step_data_long$level)) +
      labs(title = "Step Plot: Index vs. Focus Level", x = focus_label, y = "Index (Relative Scale)")

    if (panels) {
      # Faceted plot
      p <- p +
        # geom_line(aes(group = step_name), color="grey") + # Base lines (optional)
        geom_line(aes(group = step_name), color = "black") + # Current step black
        geom_point(aes(color = step_name), na.rm = TRUE) + # Color points by step, handle NAs
        facet_wrap(~step_name, ncol = 1, strip.position = "top") +
        scale_color_brewer("Step", palette = "Paired") + # Use color palette
        theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside")
    } else { # Overlaid plot
      p <- p +
        geom_line(aes(color = step_name, linetype = step_name), na.rm = TRUE) +
        geom_point(aes(color = step_name, shape = step_name), na.rm = TRUE) +
        scale_color_brewer("Step", palette = "Paired") +
        scale_linetype_discrete("Step") +
        scale_shape_discrete("Step") +
        theme(legend.position = "right") # Move legend
    }
    return(p)
  }

  # Influence Plot
  else if (type == "influ") {
    influ_data_long <- prepare_influ_data(obj$influences, focus_var)
    if (is.null(influ_data_long)) stop("No influence columns found in obj$influences.")

    ylab <- "Influence (Partial Effect on Link Scale)"

    p <- ggplot(influ_data_long, aes(x = level_num, y = influence_value)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      scale_x_continuous(breaks = unique(influ_data_long$level_num), labels = levels(influ_data_long$level)) +
      labs(title = "Influence Plot: Term Effect vs. Focus Level", x = focus_label, y = ylab)

    if (panels) {
      p <- p +
        geom_line(color = "black", na.rm = TRUE) +
        geom_point(aes(color = term_label), na.rm = TRUE) + # Color points by term
        facet_wrap(~term_label, ncol = 1, scales = "free_y", strip.position = "top") +
        scale_color_brewer("Term", palette = "Set3") + # Use color palette
        theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside")
    } else { # Overlaid plot
      p <- p +
        geom_line(aes(color = term_label, linetype = term_label), na.rm = TRUE) +
        geom_point(aes(color = term_label, shape = term_label), na.rm = TRUE) +
        scale_color_brewer("Term", palette = "Set3") +
        scale_linetype_discrete("Term") +
        scale_shape_discrete("Term") +
        theme(legend.position = "right") # Move legend
    }
    return(p)
  }

  # CDI Plot (Experimental for GAMs)
  else if (type == "cdi") {
    if (is.null(term)) stop("Argument 'term' must be provided for type = 'cdi'.")
    if (!term %in% names(obj$influences)) {
      # Try matching term labels from summary which might have '+'
      clean_summary_terms <- sub("^\\+ ", "", obj$summary$term)
      clean_summary_terms <- clean_summary_terms[clean_summary_terms != "intercept"]
      if (term %in% clean_summary_terms) {
        # Use the name found in influences if available
        matching_influ_term <- names(obj$influences)[names(obj$influences) == term]
        if (length(matching_influ_term) == 0) {
          # Check if term is part of influences names (e.g. s(x) vs s(x) in summary)
          possible_matches <- names(obj$influences)[startsWith(names(obj$influences), term)]
          if (length(possible_matches) == 1) {
            term <- possible_matches
            warning("Using term '", term, "' found in influences matching CDI term '", term, "'.", call. = FALSE)
          } else {
            stop("Term '", term, "' (or variation) not found in calculated influences. Available: ", paste(names(obj$influences)[-1], collapse = ", "))
          }
        }
      } else {
        stop("Term '", term, "' not found in model summary terms used for influences. Available: ", paste(names(obj$influences)[-1], collapse = ", "))
      }
    }
    if (is.null(obj$preds)) stop("Predictions (obj$preds) required for CDI plot are missing.")
    if (!term %in% names(obj$preds)) stop("Term '", term, "' not found in predicted term effects (obj$preds).")


    message("CDI plot for GAMs is experimental. 'Coefficient' plot shows partial effect on link scale.")

    cdi_data <- prepare_cdi_data(obj, term)
    term_label <- labels[[cdi_data$term_var]] %||% cdi_data$term_var
    focus_label <- labels[[focus_var]] %||% focus_var

    # 1. Coefficient Plot (Partial Effect vs Term Variable)
    p_coeff <- ggplot(cdi_data$coeffs, aes(x = term_var_num, y = term_effect)) +
      geom_point(size = 2, color = RColorBrewer::brewer.pal(3, "Set1")[1], na.rm = TRUE) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      labs(y = paste("Partial Effect"), x = NULL, title = paste("CDI Plot:", term)) + # Add title here
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(t = 5, r = 5, b = 0, l = 5)) # Adjust margins
    if (cdi_data$is_numeric) {
      # Optionally add a smooth line for numeric terms
      p_coeff <- p_coeff + geom_smooth(method = "loess", se = FALSE, color = "darkgrey", linetype = "dotted", na.rm = TRUE)
    } else { # Factor term
      p_coeff <- p_coeff + scale_x_continuous(breaks = cdi_data$coeffs$term_var_num, labels = cdi_data$coeffs$term_var_labels) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate labels if factor
    }

    # 2. Distribution Plot (Focus Level vs Term Variable, point size = proportion)
    p_distr <- ggplot(cdi_data$distr, aes(x = term_int, y = focus_int)) +
      geom_point(aes(size = prop), shape = 16, show.legend = TRUE, color = RColorBrewer::brewer.pal(3, "Set1")[2], na.rm = TRUE) +
      scale_size_continuous(range = c(1, 8), name = "Proportion") + # Adjust size range
      scale_x_continuous(breaks = unique(cdi_data$distr$term_int), labels = levels(cdi_data$distr$term_labels)) +
      scale_y_continuous(breaks = unique(cdi_data$distr$focus_int), labels = levels(cdi_data$distr$focus_labels)) +
      labs(x = term_label, y = focus_label) +
      theme(
        panel.grid.major = element_line(color = "grey90"),
        legend.position = "bottom",
        axis.text.x = if (!cdi_data$is_numeric) element_text(angle = 45, hjust = 1) else element_text(), # Rotate labels if factor
        plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
      ) # Adjust margins


    # 3. Influence Plot (Term Influence vs Focus Level)
    p_influ <- ggplot(cdi_data$influ, aes(x = influence_value, y = level_num)) +
      geom_point(size = 2, color = RColorBrewer::brewer.pal(3, "Set1")[3], na.rm = TRUE) +
      geom_path(color = RColorBrewer::brewer.pal(3, "Set1")[3], na.rm = TRUE) + # Connect points
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
      scale_y_continuous(breaks = cdi_data$influ$level_num, labels = cdi_data$influ$level) +
      labs(x = "Influence", y = NULL) + # Simpler label
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 5, l = 0)
      ) # Adjust margins

    # Combine using patchwork
    layout <- "
    AA#
    BBC
    BBC
    "
    combined_plot <- p_coeff + p_distr + p_influ + plot_layout(design = layout) & theme(legend.position = "bottom")
    # Add overall title using patchwork annotation
    # combined_plot <- combined_plot + plot_annotation(title = paste("CDI Plot:", term))

    return(combined_plot)
  }

  # Step + Influence Plot
  else if (type == "step_influ") {
    # Generate individual plots (faceted)
    p_step <- plot.influence_gam(x, type = "step", panels = TRUE, labels = labels, base_size = base_size, ...) +
      theme(legend.position = "none", plot.margin = margin(5, 0, 5, 5)) + labs(title = NULL, x = NULL)
    p_influ <- plot.influence_gam(x, type = "influ", panels = TRUE, labels = labels, base_size = base_size, ...) +
      theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(5, 5, 5, 0)) + labs(title = NULL, y = NULL, x = NULL)

    # Combine side-by-side
    combined_plot <- p_step + p_influ + plot_layout(ncol = 2)
    # Add overall title and shared x-axis label
    combined_plot <- combined_plot +
      plot_annotation(
        title = "Step and Influence Plots",
        caption = paste("Focus Variable:", focus_label),
        theme = theme(plot.title = element_text(hjust = 0.5))
      ) &
      theme(plot.margin = margin(1, 1, 1, 1)) # Adjust margins for subplots slightly


    return(combined_plot)
  }

  # Interaction Plot
  else if (type == "interaction") {
    if (is.null(term)) stop("Argument 'term' must be provided for type = 'interaction'.")

    interaction_data <- prepare_interaction_data(obj, term)
    plot_data <- interaction_data$data
    var1 <- interaction_data$var1
    var2 <- interaction_data$var2
    term_name <- interaction_data$term

    var1_label <- labels[[var1]] %||% var1
    var2_label <- labels[[var2]] %||% var2
    response_label <- "Interaction Effect (Link Scale)"

    p <- ggplot(plot_data, aes(
      x = x_factor, y = response_vals,
      color = trace_factor, shape = trace_factor, linetype = trace_factor, group = trace_factor
    )) +
      geom_point(size = 2.5, na.rm = TRUE) +
      geom_line(na.rm = TRUE) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      scale_color_brewer(var2_label, palette = "Set1") + # Use color palette
      scale_shape_discrete(var2_label) +
      scale_linetype_discrete(var2_label) +
      labs(
        title = paste("Interaction Plot:", term_name),
        x = var1_label, y = response_label
      ) +
      theme(legend.position = "top")
    return(p)
  }

  # Should not be reached if type is matched, but return NULL invisibly if it is
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Plotting Helper Functions (Consider placing in a separate plot-helpers.R file)
# -----------------------------------------------------------------------------

# Helper function for NULL or value
`%||%` <- function(a, b) if (is.null(a)) b else a

# Helper function to generate prediction grid for interactions
create_interaction_grid <- function(model, data, var1, var2) {
  if (!var1 %in% names(data) || !var2 %in% names(data)) {
    stop("Interaction variables (", var1, ", ", var2, ") not found in data.")
  }

  # Create grid using unique values or factor levels
  var1_vals <- if (is.factor(data[[var1]])) levels(data[[var1]]) else sort(unique(data[[var1]]))
  var2_vals <- if (is.factor(data[[var2]])) levels(data[[var2]]) else sort(unique(data[[var2]]))

  grid <- expand.grid(list(var1_vals, var2_vals), stringsAsFactors = FALSE) # Avoid auto factor conversion initially
  names(grid) <- c(var1, var2)

  # Add reference values for other predictors
  model_vars <- all.vars(stats::formula(model)[-2]) # All predictors in formula RHS
  other_vars <- setdiff(model_vars, c(var1, var2))

  for (v in other_vars) {
    if (!v %in% names(data)) next # Skip if var not in data (e.g., part of smooth spec like 'k')
    if (is.numeric(data[[v]])) {
      grid[[v]] <- mean(data[[v]], na.rm = TRUE)
    } else if (is.factor(data[[v]])) {
      # Use the first level as the reference level
      grid[[v]] <- factor(levels(data[[v]])[1], levels = levels(data[[v]]))
    } else if (is.character(data[[v]])) {
      # If character, convert to factor using levels from original data
      grid[[v]] <- factor(data[[v]][1], levels = unique(data[[v]]))
    } else {
      grid[[v]] <- data[[v]][1] # Fallback for other types (logical, etc.)
    }
  }

  # Explicitly set factor types for interaction vars based on original data
  if (is.factor(data[[var1]])) grid[[var1]] <- factor(grid[[var1]], levels = levels(data[[var1]]))
  if (is.factor(data[[var2]])) grid[[var2]] <- factor(grid[[var2]], levels = levels(data[[var2]]))

  return(grid)
}


# Helper function for step plot data preparation
prepare_step_data <- function(indices_df, focus_var) {
  # Identify step columns (heuristic: not level, unstan, stan, stanLower, stanUpper)
  step_cols <- setdiff(names(indices_df), c("level", "unstan", "stan", "stanLower", "stanUpper"))
  if (length(step_cols) == 0) {
    return(NULL)
  }

  step_data <- indices_df[, c("level", step_cols), drop = FALSE]
  step_data_long <- tryCatch(
    tidyr::pivot_longer(step_data,
      cols = all_of(step_cols), # Use all_of for safety
      names_to = "step_name",
      values_to = "index_value"
    ),
    error = function(e) {
      warning("tidyr::pivot_longer failed for step data: ", e$message)
      NULL
    }
  )
  if (is.null(step_data_long)) {
    return(NULL)
  }

  # Ensure step_name is ordered correctly based on original column order
  step_data_long$step_name <- factor(step_data_long$step_name, levels = step_cols)
  # Ensure level is ordered correctly
  step_data_long$level <- factor(step_data_long$level, levels = indices_df$level)
  step_data_long$level_num <- as.numeric(step_data_long$level)
  return(step_data_long)
}

# Helper function for influence plot data preparation
prepare_influ_data <- function(influences_df, focus_var) {
  # Identify influence columns (all except 'level')
  influ_cols <- setdiff(names(influences_df), "level")
  if (length(influ_cols) == 0) {
    return(NULL)
  }

  influ_data <- influences_df[, c("level", influ_cols), drop = FALSE]
  influ_data_long <- tryCatch(
    tidyr::pivot_longer(influ_data,
      cols = all_of(influ_cols), # Use all_of
      names_to = "term_label",
      values_to = "influence_value"
    ),
    error = function(e) {
      warning("tidyr::pivot_longer failed for influence data: ", e$message)
      NULL
    }
  )
  if (is.null(influ_data_long)) {
    return(NULL)
  }

  # Ensure term_label is ordered correctly
  influ_data_long$term_label <- factor(influ_data_long$term_label, levels = influ_cols)
  # Ensure level is ordered correctly
  influ_data_long$level <- factor(influ_data_long$level, levels = influences_df$level)
  influ_data_long$level_num <- as.numeric(influ_data_long$level)
  return(influ_data_long)
}


# Helper function for CDI plot data preparation
prepare_cdi_data <- function(obj, term) {
  preds <- obj$preds
  influences <- obj$influences
  data <- obj$data
  focus_var <- obj$focus_var

  # Identify the primary variable associated with the term (heuristic for axis label)
  # Handles s(x), s(x,k=...), ti(x,y), x:y, x
  term_var <- gsub("^s\\(|^te\\(|^ti\\(", "", term) # Remove smooth wrappers
  term_var <- gsub(",.*$", "", term_var) # Remove args like k=...
  term_var <- gsub("\\)$", "", term_var) # Remove trailing )
  term_var <- strsplit(term_var, ":")[[1]][1] # Take first part of interaction
  term_var <- trimws(term_var)

  if (!term_var %in% names(data)) {
    # Fallback or warning if heuristic fails
    warning("Could not robustly determine primary variable for term '", term, "'. Using term name for axis label.", call. = FALSE)
    term_var_for_labeling <- term # Use full term name if base variable not found
  } else {
    term_var_for_labeling <- term_var
  }

  # Get the actual values of the variable(s) involved in the term from original data
  # If it's a smooth term, we often want the first variable mentioned
  actual_term_var <- term_var # Start with the guess
  if (!actual_term_var %in% names(data)) {
    # If primary var not found, maybe it's implicit (e.g. factor name used directly)
    if (term %in% names(data)) {
      actual_term_var <- term
      term_var_for_labeling <- term
    } else {
      stop("Cannot find variable(s) in data corresponding to term '", term, "' for CDI plot.")
    }
  }


  plot_data <- data.frame(
    term_effect = preds[[term]],
    term_var_values = data[[actual_term_var]],
    focus_var_values = data[[focus_var]]
  )
  plot_data <- na.omit(plot_data) # Important: remove rows with NAs in relevant columns

  # Data for Coefficient plot (Partial Effect vs Term Variable)
  is_numeric_term <- is.numeric(plot_data$term_var_values)
  if (is_numeric_term) {
    # Use original numeric values for plotting, maybe aggregate for points if too many
    if (nrow(plot_data) > 500) { # Aggregate if very large dataset
      breaks <- pretty(plot_data$term_var_values, n = 20) # Create bins
      if (length(breaks) < 2) breaks <- range(plot_data$term_var_values)
      midpoints <- (utils::head(breaks, -1) + utils::tail(breaks, -1)) / 2
      plot_data$term_var_cut <- cut(plot_data$term_var_values, breaks, labels = midpoints, include.lowest = TRUE)
      coeffs_data <- stats::aggregate(term_effect ~ term_var_cut, data = plot_data, mean, na.rm = TRUE)
      coeffs_data$term_var_num <- as.numeric(as.character(coeffs_data$term_var_cut)) # Midpoints as numeric x
      coeffs_data$term_var_labels <- coeffs_data$term_var_cut # Keep labels if needed for hover?
    } else { # Use raw points if not too many
      coeffs_data <- plot_data[, c("term_var_values", "term_effect")]
      names(coeffs_data) <- c("term_var_num", "term_effect")
      coeffs_data$term_var_labels <- as.character(coeffs_data$term_var_num)
    }
  } else { # Factor term
    coeffs_data <- stats::aggregate(term_effect ~ term_var_values, data = plot_data, mean, na.rm = TRUE)
    # Ensure factor levels match original data for correct ordering
    coeffs_data$term_var_values <- factor(coeffs_data$term_var_values, levels = levels(data[[actual_term_var]]))
    coeffs_data <- coeffs_data[order(coeffs_data$term_var_values), ] # Order by factor level
    coeffs_data$term_var_num <- as.numeric(coeffs_data$term_var_values) # Numeric position
    coeffs_data$term_var_labels <- coeffs_data$term_var_values # Factor labels
  }


  # Data for Distribution plot (Focus level vs Term variable)
  distr_data <- as.data.frame(table(plot_data$term_var_values, plot_data$focus_var_values, dnn = c("term_val", "focus_val")))
  # Use original levels for factors
  if (is.factor(data[[actual_term_var]])) distr_data$term_val <- factor(distr_data$term_val, levels = levels(data[[actual_term_var]]))
  if (is.factor(data[[focus_var]])) distr_data$focus_val <- factor(distr_data$focus_val, levels = levels(data[[focus_var]]))

  distr_data <- subset(distr_data, Freq > 0)
  distr_data$term_int <- as.integer(factor(distr_data$term_val)) # Consistent numeric position
  distr_data$focus_int <- as.integer(factor(distr_data$focus_val)) # Consistent numeric position

  total_per_focus <- stats::aggregate(Freq ~ focus_val, data = distr_data, sum)
  distr_data <- merge(distr_data, total_per_focus, by = "focus_val", suffixes = c("", "_total"))
  distr_data$prop <- distr_data$Freq / distr_data$Freq_total

  # Get labels matching the integer positions
  term_labels_map <- unique(distr_data[, c("term_int", "term_val")])
  focus_labels_map <- unique(distr_data[, c("focus_int", "focus_val")])
  distr_data$term_labels <- term_labels_map$term_val[match(distr_data$term_int, term_labels_map$term_int)]
  distr_data$focus_labels <- focus_labels_map$focus_val[match(distr_data$focus_int, focus_labels_map$focus_int)]
  # Ensure factors use original levels for plotting order
  if (is.factor(distr_data$term_labels)) distr_data$term_labels <- factor(distr_data$term_labels, levels = levels(data[[actual_term_var]]))
  if (is.factor(distr_data$focus_labels)) distr_data$focus_labels <- factor(distr_data$focus_labels, levels = levels(data[[focus_var]]))


  # Data for Influence plot (Influence value vs Focus Level)
  if (!term %in% names(influences)) {
    # If the exact term name (e.g., 's(x)') isn't in influences, check if a related one is
    possible_matches <- names(influences)[grepl(term_var_for_labeling, names(influences), fixed = TRUE)]
    if (length(possible_matches) == 1) {
      term_in_influences <- possible_matches
      warning("Using influence term '", term_in_influences, "' for CDI plot based on base variable '", term_var_for_labeling, "'.", call. = FALSE)
    } else {
      stop("Cannot find unique influence term related to '", term, "' for CDI plot. Found: ", paste(possible_matches, collapse = ", "))
    }
  } else {
    term_in_influences <- term
  }

  influ_plot_data <- data.frame(
    level = influences$level,
    influence_value = influences[[term_in_influences]]
  )
  # Ensure level is factor with correct order
  influ_plot_data$level <- factor(influ_plot_data$level, levels = levels(data[[focus_var]]))
  influ_plot_data <- influ_plot_data[order(influ_plot_data$level), ] # Order by level
  influ_plot_data$level_num <- as.numeric(influ_plot_data$level)


  return(list(
    coeffs = coeffs_data,
    distr = distr_data,
    influ = influ_plot_data,
    term_var = term_var_for_labeling, # Use potentially simplified name for label
    is_numeric = is_numeric_term
  ))
}


# Helper function for Interaction plot data preparation
prepare_interaction_data <- function(obj, term) {
  interaction_vars <- strsplit(term, ":")[[1]]
  if (length(interaction_vars) != 2) {
    stop("Interaction plot currently only supports two-way interaction terms specified as 'var1:var2'.")
  }
  var1 <- trimws(interaction_vars[1])
  var2 <- trimws(interaction_vars[2])

  if (!var1 %in% names(obj$data) || !var2 %in% names(obj$data)) {
    stop("Variables '", var1, "' or '", var2, "' not found in the model data.")
  }
  if (!is.factor(obj$data[[var1]]) || !is.factor(obj$data[[var2]])) {
    warning("Interaction plot works best with two factor variables. Results may be unexpected.", call. = FALSE)
    # Allow non-factors, but plotting might be odd
  }

  pred_grid <- create_interaction_grid(obj$model, obj$data, var1, var2)

  all_term_preds <- tryCatch(stats::predict(obj$model, newdata = pred_grid, type = "terms", se.fit = FALSE),
    error = function(e) {
      stop("Failed to predict terms for interaction plot: ", e$message)
    }
  )

  # Check for term name variations (e.g., var1:var2 vs var2:var1)
  term_to_use <- NULL
  if (term %in% colnames(all_term_preds)) {
    term_to_use <- term
  } else {
    term_reversed <- paste(var2, var1, sep = ":")
    if (term_reversed %in% colnames(all_term_preds)) {
      term_to_use <- term_reversed
      message("Using term '", term_to_use, "' found in predictions for interaction plot.")
    } else {
      # Handle potential smooth interaction names (e.g., ti(var1,var2))
      smooth_term_pattern <- paste0("(ti|te|t2)\\(.*", var1, ".*,", ".*", var2, ".*\\)")
      smooth_term_pattern_rev <- paste0("(ti|te|t2)\\(.*", var2, ".*,", ".*", var1, ".*\\)")
      matching_smooth <- grep(smooth_term_pattern, colnames(all_term_preds), value = TRUE)
      matching_smooth_rev <- grep(smooth_term_pattern_rev, colnames(all_term_preds), value = TRUE)

      if (length(matching_smooth) == 1) {
        term_to_use <- matching_smooth
      } else if (length(matching_smooth_rev) == 1) term_to_use <- matching_smooth_rev

      if (is.null(term_to_use)) {
        stop("Specified interaction term '", term, "' (or variations) not found in predicted terms: ", paste(colnames(all_term_preds), collapse = ", "))
      } else {
        message("Using smooth term '", term_to_use, "' found in predictions for interaction plot.")
      }
    }
  }

  pred_grid$interaction_effect <- all_term_preds[, term_to_use]

  # Decide which variable goes on x-axis vs trace factor (var1 on x, var2 as trace)
  pred_grid$x_factor <- pred_grid[[var1]]
  pred_grid$trace_factor <- pred_grid[[var2]]
  pred_grid$response_vals <- pred_grid$interaction_effect

  # Ensure factors have correct levels if they were character originally
  if (is.factor(obj$data[[var1]])) pred_grid$x_factor <- factor(pred_grid$x_factor, levels = levels(obj$data[[var1]]))
  if (is.factor(obj$data[[var2]])) pred_grid$trace_factor <- factor(pred_grid$trace_factor, levels = levels(obj$data[[var2]]))


  return(list(data = pred_grid, var1 = var1, var2 = var2, term = term_to_use))
}
