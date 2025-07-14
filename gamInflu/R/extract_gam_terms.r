#' Extract GAM Terms Information
#'
#' Internal function to extract information about smooth and parametric terms
#'
#' @param model A fitted GAM model
#' @return List with term information
#' @keywords internal
#' 
extract_gam_terms <- function(model) {
  # Get smooth terms
  smooth_terms <- character(0)
  smooth_info <- list()

  if (length(model$smooth) > 0) {
    for (i in seq_along(model$smooth)) {
      smooth <- model$smooth[[i]]
      term_label <- smooth$label
      smooth_terms <- c(smooth_terms, term_label)

      smooth_info[[term_label]] <- list(
        label = term_label,
        dim = smooth$dim,
        by = smooth$by,
        bs = smooth$bs.dim,
        edf = model$edf[i],
        variables = smooth$term,
        by_variable = if (smooth$by != "NA") smooth$by else NULL
      )
    }
  }

  # Get parametric terms
  param_terms <- character(0)
  if (length(model$coefficients) > length(unlist(lapply(model$smooth, function(x) x$label)))) {
    param_labels <- names(model$coefficients)
    smooth_labels <- unlist(lapply(model$smooth, function(x) x$label))
    param_terms <- setdiff(param_labels, smooth_labels)
  }

  return(list(
    smooth_terms = smooth_terms,
    parametric_terms = param_terms,
    smooth_info = smooth_info,
    n_smooth = length(smooth_terms),
    n_parametric = length(param_terms)
  ))
}
