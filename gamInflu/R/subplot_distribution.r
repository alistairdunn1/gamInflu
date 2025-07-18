#' @title Plot Data Distribution for a Model Term
#' @description Plots the data distribution for a specified model term and focus variable.
#' @param obj A `gam_influence` object.
#' @param term The character name of the model term to plot (e.g., "year" or "s(temp)").
#' @return A ggplot object showing the data distribution.
#' @export
plot_term_distribution <- function(obj, term) {
  # --- Setup ---
  if (is.numeric(term) && length(term) == 1 && term == as.integer(term)) {
    all_terms <- get_terms(obj, full = TRUE)
    if (term > length(all_terms) || term < 1) {
      stop("Term index out of bounds. There are ", length(all_terms), " valid terms:\n  ", paste(all_terms, collapse = "\n  "), call. = FALSE)
    }
    term <- all_terms[term]
  }

  # plot
  obj_data <- obj$data
  focus_var <- obj$focus
  if (term == focus_var) {
    stop("The distribution plot cannot be generated for the focus term itself.", call. = FALSE)
  }
  p_dist <- subplot_distribution(obj_data, term, focus_var)
  return(p_dist)
}

#' @title Subplot for Data Distribution
#' @description Internal function to plot the data distribution for a term in CDI plots.
#' @param p_dist_data Data frame with distribution counts.
#' @param focus_var The focus variable name.
#' @param term_var The term variable name.
#' @noRd
subplot_distribution <- function(obj_data, term, focus_var) {
  term_vars <- all.vars(rlang::parse_expr(term))
  call_obj <- rlang::parse_expr(term)
  call_list <- as.list(call_obj)

  if (is.numeric(obj_data[[term_vars[1]]])) {
    # Cut numeric variable into 10 factor levels and label by midpoint
    cuts <- cut(obj_data[[term_vars[1]]], breaks = 10, include.lowest = TRUE)
    # Get midpoints for each interval
    cut_levels <- levels(cuts)
    cut_midpoints <- sapply(strsplit(gsub("\\[|\\]|\\(|\\)", "", cut_levels), ","), function(x) {
      mean(as.numeric(x))
    })
    # Relabel factor levels with midpoints
    levels(cuts) <- round(cut_midpoints, 2)
    obj_data[[term_vars[1]]] <- cuts
  }
  p_dist_data <- obj_data %>%
    dplyr::group_by(.data[[focus_var]], .data[[term_vars[1]]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      n = ((n * 8) / max(n, na.rm = TRUE)),
      n = if (length(unique(n)) == 1) 1 else n
    )


  p_dist <- ggplot(p_dist_data, aes(x = !!rlang::sym(focus_var), y = !!rlang::sym(term_vars[1]))) +
    geom_point(colour = "royalblue", size = p_dist_data$n, alpha = 0.5) +
    coord_flip() +
    labs(y = term_vars[1], x = focus_var)

  return(p_dist)
}
