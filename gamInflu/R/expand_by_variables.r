#' Expand By-Variables into Separate Terms
#'
#' @param model GAM model
#' @return List with expanded term information
#' @keywords internal
#'
expand_by_variables <- function(model) {
  expanded_terms <- list()
  used_labels <- character(0) # Track used labels to avoid conflicts

  # Get all terms in the model formula
  all_terms <- attr(terms(formula(model)), "term.labels")
  smooth_labels <- sapply(model$smooth, function(s) s$label)

  for (i in seq_along(model$smooth)) {
    smooth <- model$smooth[[i]]

    if (smooth$by != "NA") {
      # This is a by-variable smooth
      by_var <- smooth$by
      by_levels <- levels(model$model[[by_var]])

      # Create separate terms for each by-level
      for (level in by_levels) {
        term_label <- level # Just use the by-level name

        # Check for naming conflicts and resolve them
        if (term_label %in% used_labels) {
          # Add original smooth info to distinguish
          term_label <- paste0(gsub("[^a-zA-Z0-9]", "_", smooth$label), "_", level)
          warning("Naming conflict resolved: using '", term_label, "' for by-level '", level, "'")
        }

        used_labels <- c(used_labels, term_label)

        expanded_terms[[term_label]] <- list(
          original_label = smooth$label,
          by_variable = by_var,
          by_level = level,
          smooth_index = i,
          variables = smooth$term,
          bs_dim = smooth$bs.dim,
          type = "by_smooth"
        )
      }
    } else {
      # Regular smooth term
      term_label <- smooth$label

      # Check for conflicts with by-variable labels
      if (term_label %in% used_labels) {
        term_label <- paste0("smooth_", term_label)
        warning("Naming conflict resolved: using '", term_label, "' for smooth term")
      }

      used_labels <- c(used_labels, term_label)

      expanded_terms[[term_label]] <- list(
        original_label = smooth$label,
        smooth_index = i,
        variables = smooth$term,
        bs_dim = smooth$bs.dim,
        type = "smooth"
      )
    }
  }

  # Add regular (parametric) terms: those in the formula but not in smooths
  parametric_terms <- setdiff(all_terms, smooth_labels)
  for (term in parametric_terms) {
    # Check for conflicts
    term_label <- term
    if (term_label %in% used_labels) {
      term_label <- paste0("param_", term_label)
      warning("Naming conflict resolved: using '", term_label, "' for parametric term")
    }
    used_labels <- c(used_labels, term_label)

    expanded_terms[[term_label]] <- list(
      original_label = term,
      variables = term,
      type = "parametric"
    )
  }

  return(expanded_terms)
}
