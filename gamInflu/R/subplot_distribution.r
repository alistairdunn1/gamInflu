#' @title Plot Data Distribution for a Model Term
#' @description Plots the data distribution for a specified model term and focus variable.
#' @param obj A `gam_influence` object.
#' @param term The character name of the model term to plot (e.g., "year" or "s(temp)").
#' @param by Logical indicating whether to plot by-variable interactions (default is TRUE).
#'   If FALSE, it will plot the main variable only.
#' @return A ggplot object showing the data distribution.
#' @export
plot_term_distribution <- function(obj, term, by = TRUE) {
  message("Plotting data distribution for term: ", term)
  # --- Setup ---
  if (is.numeric(term) && length(term) == 1 && term == as.integer(term)) {
    all_terms <- get_terms(obj, full = TRUE)
    if (term > length(all_terms) || term < 1) {
      stop("Term index out of bounds. There are ", length(all_terms), " valid terms:\n  ", paste(all_terms, collapse = "\n  "))
    }
    term <- all_terms[term]
  }

  # plot
  focus_var <- obj$focus
  if (term == focus_var) {
    stop("The distribution plot cannot be generated for the focus term itself.")
  }
  p_dist <- subplot_distribution(obj, term, focus_var, by = by)
  return(p_dist)
}

#' @title Subplot for Data Distribution
#' @description Internal function to plot the data distribution for a term in CDI plots.
#' @param obj A `gam_influence` object.
#' @param term The character name of the model term to plot (e.g., "year" or "s(temp)").
#' @param focus_var The focus variable name.
#' @param by Logical indicating whether to plot by-variable interactions (default is TRUE).
#' @return A ggplot object showing the data distribution.
#' @noRd
subplot_distribution <- function(obj, term, focus_var, by = TRUE) {
  term_vars <- all.vars(rlang::parse_expr(term))

  # For interactions, we need to handle multiple variables
  # First, identify which variables are actually in the model terms
  model_terms <- get_terms(obj, full = FALSE)
  term_vars_in_model <- term_vars[term_vars %in% model_terms]

  if (length(term_vars_in_model) == 0) {
    stop(paste0(
      "None of the variables in the supplied term ('", term, "') match the model terms.\n",
      "Valid terms are: ", paste(setdiff(model_terms, obj$focus), collapse = ", ")
    ))
  }

  # For multiple variables (interactions), use the first variable in the model for display
  # but consider this could be enhanced to show interaction structure
  main_var <- term_vars_in_model[1]

  # Validate that focus_var and main_var are valid character strings
  if (!is.character(focus_var) || length(focus_var) != 1 || focus_var == "" || is.na(focus_var)) {
    stop("focus_var must be a valid character string")
  }
  if (!is.character(main_var) || length(main_var) != 1 || main_var == "" || is.na(main_var)) {
    stop("main_var must be a valid character string")
  }

  # Check if this is a by-variable interaction (e.g., s(year, by=area))
  is_by_interaction <- grepl("by\\s*=", term)

  # Handle binning for numeric variables
  if (is.numeric(obj$data[[main_var]])) {
    n <- unique(obj$data[[main_var]])
    if (length(n) > 15) {
      breaks <- 10
      cuts <- cut(obj$data[[main_var]], breaks = breaks, include.lowest = TRUE)
      # Get midpoints for each interval
      cut_levels <- levels(cuts)
      cut_midpoints <- sapply(strsplit(gsub("\\[|\\]|\\(|\\)", "", cut_levels), ","), function(x) {
        mean(as.numeric(x))
      })
      # Relabel factor levels with midpoints
      levels(cuts) <- round(cut_midpoints, 2)
      obj$data[[main_var]] <- cuts
    } else {
      # For 15 or fewer unique values, use the actual values as factor levels
      obj$data[[main_var]] <- factor(obj$data[[main_var]], levels = sort(n))
    }
  }

  # Handle different types of interactions
  if (length(term_vars_in_model) > 1 && !is_by_interaction) {
    # For tensor product or other multi-variable interactions (e.g., te(year, area))
    second_var <- term_vars_in_model[2]

    # Validate second_var
    if (!is.character(second_var) || length(second_var) != 1 || second_var == "" || is.na(second_var)) {
      stop("second_var must be a valid character string")
    }

    p_dist_data <- obj$data %>%
      dplyr::group_by(.data[[focus_var]], .data[[main_var]], .data[[second_var]]) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        n = ((n * 8) / max(n, na.rm = TRUE)),
        n = if (length(unique(n)) == 1) 1 else n
      )

    p_dist <- ggplot(p_dist_data, aes(x = !!rlang::sym(focus_var), y = !!rlang::sym(main_var))) +
      geom_point(aes(size = n), colour = "royalblue", alpha = 0.5) +
      ggplot2::facet_wrap(~ !!rlang::sym(second_var), scales = "free") +
      ggplot2::coord_flip() +
      labs(y = paste(main_var, "x", second_var), x = focus_var) +
      ggplot2::scale_size_identity() +
      ggplot2::theme(strip.text = ggplot2::element_text(size = ggplot2::rel(0.8)))
  } else if (is_by_interaction) {
    if (by == TRUE) {
      # For by-variable interactions (e.g., s(year, by=area))
      if (!is.character(term) || length(term) != 1 || term == "" || is.na(term)) {
        stop("term must be a valid character string")
      }
      # For by-variable interactions (e.g., s(year, by=area))
      by_var <- sub(".*by\\s*=\\s*([^,\\)]+).*", "\\1", term)
      # Clean the by_var string by removing any whitespace and quotes
      by_var <- trimws(gsub("[\"']", "", by_var))

      # Ensure by_var is a valid character string
      if (!is.character(by_var) || length(by_var) != 1 || by_var == "" || is.na(by_var)) {
        warning("Could not extract valid by-variable from term: ", term)
        by_var <- NULL
      }

      # Additional validation: check if by_var is a valid R identifier
      if (!is.null(by_var) && !grepl("^[a-zA-Z][a-zA-Z0-9_]*$", by_var)) {
        warning("by_var '", by_var, "' is not a valid R identifier")
        by_var <- NULL
      }

      p_dist_data <- obj$data %>%
        dplyr::group_by(.data[[focus_var]], .data[[main_var]], .data[[by_var]]) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          n = ((n * 8) / max(n, na.rm = TRUE)),
          n = if (length(unique(n)) == 1) 1 else n
        )

      # Create the plot with error handling around the facet_wrap call
      p_dist <- ggplot(p_dist_data, aes(x = !!rlang::sym(focus_var), y = !!rlang::sym(main_var))) +
        geom_point(aes(size = n), colour = "royalblue", alpha = 0.5) +
        ggplot2::coord_flip() +
        labs(y = paste(main_var, "by", by_var), x = focus_var) +
        ggplot2::scale_size_identity() +
        ggplot2::theme(strip.text = ggplot2::element_text(size = ggplot2::rel(0.8)))

      # Create a formula for faceting using the variable name directly
      facet_formula <- as.formula(paste("~", by_var))
      p_dist <- p_dist + ggplot2::facet_wrap(facet_formula, scales = "free")
    } else {
      # Fallback to single variable if by-variable not found
      p_dist_data <- obj$data %>%
        dplyr::group_by(.data[[focus_var]], .data[[main_var]]) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          n = ((n * 8) / max(n, na.rm = TRUE)),
          n = if (length(unique(n)) == 1) 1 else n
        )

      p_dist <- ggplot(p_dist_data, aes(x = !!rlang::sym(focus_var), y = !!rlang::sym(main_var))) +
        geom_point(aes(size = n), colour = "royalblue", alpha = 0.5) +
        ggplot2::coord_flip() +
        labs(y = main_var, x = focus_var) +
        ggplot2::scale_size_identity()
    }
  } else {
    # Single variable case
    p_dist_data <- obj$data %>%
      dplyr::group_by(.data[[focus_var]], .data[[main_var]]) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        n = ((n * 8) / max(n, na.rm = TRUE)),
        n = if (length(unique(n)) == 1) 1 else n
      )

    p_dist <- ggplot(p_dist_data, aes(x = !!rlang::sym(focus_var), y = !!rlang::sym(main_var))) +
      geom_point(aes(size = n), colour = "royalblue", alpha = 0.5) +
      ggplot2::coord_flip() +
      labs(y = main_var, x = focus_var) +
      ggplot2::scale_size_identity()
  }
  return(p_dist)
}
