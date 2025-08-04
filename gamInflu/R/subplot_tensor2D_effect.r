#' @title Subplot for 2D Tensor Smooth Effects
#' @description Internal function to plot predicted effects for 2D tensor smooth terms (e.g., te(x, y)).
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param term_vars Variables in the term.
#' @param cdi Logical; if TRUE, use CDI plot style.
#' @noRd
subplot_tensor2d_effect <- function(obj, t, term_vars, cdi = FALSE) {
  message("Plotting 2D tensor effects for term: ", t)
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se
  islog <- isTRUE(obj$islog)

  zvar <- names(preds_df)[
    vapply(names(preds_df), function(nm) {
      grepl(term_vars[1], nm) && grepl(term_vars[2], nm)
    }, logical(1))
  ]

  # Check if we found a matching column
  if (length(zvar) == 0) {
    stop(
      "Could not find prediction column for 2D tensor term: ", t,
      "\nAvailable columns: ", paste(names(preds_df), collapse = ", ")
    )
  }

  # Use the first matching column if multiple found
  if (length(zvar) > 1) {
    message("Multiple matching columns found for tensor term, using: ", zvar[1])
    zvar <- zvar[1]
  }

  preds_df$zvar <- preds_df[[zvar]]

  if (islog) {
    preds_df$zvar <- exp(preds_df$zvar)
  }
  # Cut x and y into 30 groups each and label by midpoint
  x_cut <- cut(obj$data[[term_vars[1]]], breaks = 50, include.lowest = TRUE)
  y_cut <- cut(obj$data[[term_vars[2]]], breaks = 50, include.lowest = TRUE)

  # Relabel factor levels with midpoints for x
  x_levels <- levels(x_cut)
  x_midpoints <- sapply(strsplit(gsub("\\[|\\]|\\(|\\)", "", x_levels), ","), function(x) mean(as.numeric(x)))
  levels(x_cut) <- round(x_midpoints, 2)

  # Relabel factor levels with midpoints for y
  y_levels <- levels(y_cut)
  y_midpoints <- sapply(strsplit(gsub("\\[|\\]|\\(|\\)", "", y_levels), ","), function(x) mean(as.numeric(x)))
  levels(y_cut) <- round(y_midpoints, 2)

  # Aggregate mean z for each (x, y) group
  df <- data.frame(
    x = x_cut,
    y = y_cut,
    z = preds_df$zvar
  ) %>%
    dplyr::group_by(x, y) %>%
    dplyr::summarise(z = mean(z, na.rm = TRUE), .groups = "drop")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = z)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(begin = 0, end = 1, option = "D", na.value = "grey90", limits = if (islog) c(0, NA_real_) else NULL) +
    ggplot2::labs(x = term_vars[1], y = term_vars[2], fill = "Effect") +
    ggplot2::theme(panel.grid = ggplot2::element_blank())


  if (cdi) {
    p <- p +
      ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.height = ggplot2::unit(0.5, "lines"),
        legend.key.width = ggplot2::unit(1.5, "lines"),
        legend.text = ggplot2::element_text(size = ggplot2::rel(0.75))
      )
  }
  return(p)
}
