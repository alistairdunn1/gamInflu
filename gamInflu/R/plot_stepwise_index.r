#' @title Stepwise Index Plot
#' @description Creates a step plot showing how the index for the focus term changes as each model
#' term is added sequentially. This visualization helps understand the contribution of each term
#' to the final standardized index and works with all supported GLM families.
#' @param obj A `gam_influence` object containing calculated indices from `calculate_influence()`.
#' @param show_previous Logical; if TRUE, shows previous steps on each panel in color with a legend,
#'   allowing you to see how the index evolves as terms are added. If FALSE (default), shows only
#'   the current step for each panel.
#' @return A ggplot object with stepwise index plots for each term added. The plot shows:
#'   - **Individual panels**: One for each term added to the model in sequence
#'   - **Index evolution**: How the focus term's index changes with each model term
#'   - **Reference line**: Horizontal dashed line at y=1 for relative comparison
#'   - **Color coding**: (when show_previous=TRUE) Different colors for each step in the progression
#' @details
#' The stepwise index plot visualizes the model building process by showing how the focus term's
#' index changes as each term is sequentially added to the model. This helps identify:
#' - Which terms have the largest impact on the focus index
#' - Whether the index stabilizes or continues to change with additional terms
#' - The cumulative effect of model complexity on the final index
#'
#' **Family Support:**
#' Works with all supported GLM families (Gaussian, binomial, gamma, Poisson) and uses
#' family-appropriate index calculations throughout the stepwise process.
#'
#' **Panel Labels:**
#' The function uses step labels from the model when available, otherwise falls back to
#' term names. Labels are presented in the order terms appear in the original model formula.
#' @examples
#' \dontrun{
#' # Basic stepwise plot
#' gi <- gam_influence(your_model, focus = "year")
#' gi <- calculate_influence(gi)
#' plot_stepwise_index(gi)
#'
#' # Show progression with previous steps
#' plot_stepwise_index(gi, show_previous = TRUE)
#'
#' # Through generic plot method
#' plot(gi, type = "step", show_previous = TRUE)
#' }
#' @importFrom dplyr all_of
#' @importFrom ggplot2 ggplot aes geom_hline geom_line geom_point facet_wrap labs theme element_text ylim
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @export
plot_stepwise_index <- function(obj, show_previous = FALSE) {
  # Select only the columns representing step-wise indices
  step_cols <- names(obj$calculated$indices)[grepl("^\\+", names(obj$calculated$indices))]
  df_long <- tidyr::pivot_longer(obj$calculated$indices,
    cols = all_of(step_cols), names_to = "term", values_to = "index"
  )
  df_long$term <- factor(df_long$term, levels = unique(df_long$term)) # Preserve order

  # Convert level to numeric if possible
  if (is.factor(df_long$level) && all(!is.na(as.numeric(as.character(levels(df_long$level)))))) {
    df_long$level <- as.numeric(as.character(df_long$level))
  }

  # Use step_labels if present for facet labels, and preserve GAM model order
  label_map <- NULL
  if (!is.null(obj$calculated$step_labels)) {
    label_map <- obj$calculated$step_labels
    # Use the order of step_cols as in the model call
    ordered_labels <- label_map[step_cols]
    df_long$label <- factor(label_map[as.character(df_long$term)], levels = ordered_labels)
  } else {
    df_long$label <- factor(as.character(df_long$term), levels = step_cols)
  }

  if (show_previous) {
    term_levels <- step_cols # Use model call order for panels
    expanded <- do.call(rbind, lapply(seq_along(term_levels), function(i) {
      prev_terms <- term_levels[1:i]
      subset(df_long, term %in% prev_terms)[, c("level", "index", "term", "label")]
    }))
    expanded$shown_term <- factor(rep(term_levels, sapply(seq_along(term_levels), function(i) sum(df_long$term %in% term_levels[1:i]))),
      levels = term_levels
    )

    # Use label_map if available, otherwise use term names
    if (!is.null(label_map)) {
      expanded$shown_label <- factor(label_map[as.character(expanded$shown_term)], levels = label_map[term_levels])
    } else {
      expanded$shown_label <- factor(as.character(expanded$shown_term), levels = term_levels)
    }

    ggplot2::ggplot(expanded, ggplot2::aes(x = .data$level, y = .data$index, group = .data$term, colour = .data$term)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey") +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~shown_label, ncol = 1, scales = "fixed") +
      ggplot2::labs(x = obj$focus, y = "Index", colour = "Step") +
      ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0)) +
      ggplot2::ylim(0, NA)
  } else {
    ggplot2::ggplot(df_long, ggplot2::aes(x = .data$level, y = .data$index, group = 1)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey") +
      ggplot2::geom_line(color = "royalblue", alpha = 0.5) +
      ggplot2::geom_point(color = "royalblue") +
      ggplot2::facet_wrap(~label, ncol = 1, scales = "fixed") +
      ggplot2::labs(x = obj$focus, y = "Index") +
      ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0)) +
      ggplot2::ylim(0, NA)
  }
}
