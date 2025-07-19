#' @title Subplot for 2D Tensor Smooth Effects
#' @description Internal function to plot predicted effects for 2D tensor smooth terms (e.g., te(x, y)).
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param term_vars Variables in the term.
#' @param cdi Logical; if TRUE, use CDI plot style.
#' @noRd
subplot_tensor2d_effect <- function(obj, t, term_vars, cdi = FALSE) {
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se
  islog <- isTRUE(obj$islog)
  se_col <- if (!is.null(se_df)) grep(paste0("(^|\\(|\\:)", t, "($|\\)|\\:|\\d+)"), colnames(se_df), value = TRUE) else character(0)

  xvar <- term_vars[1]
  yvar <- term_vars[2]
  effect <- preds_df[[t]]
  se <- if (length(se_col) > 0) se_df[[se_col[1]]] else NA

  if (islog) effect <- exp(effect)

  df <- data.frame(x = preds_df[[xvar]], y = preds_df[[yvar]], effect = effect)

  p <- ggplot(df, aes(x = x, y = y, fill = effect)) +
    geom_raster() +
    scale_fill_viridis_c() +
    labs(x = xvar, y = yvar, fill = "Effect")

  if (cdi) {
    p <- p +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
      )
  }
  return(p)
}
