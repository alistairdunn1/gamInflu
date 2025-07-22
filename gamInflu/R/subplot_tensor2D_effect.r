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

  p <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    labs(x = term_vars[1], y = term_vars[2], fill = "Effect")

  if (cdi) {
    p <- p +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "inside",
        legend.justification = "top",
      )
  }
  return(p)
}
