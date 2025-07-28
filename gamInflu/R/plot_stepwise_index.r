#' @title Stepwise Index Plot
#' @description Creates a step plot showing how the index for the focus term changes as each model term is added sequentially.
#' @param obj A `gam_influence` object containing calculated indices.
#' @param show_previous Logical; if TRUE, shows previous steps on each panel in colour with a legend.
#' @return A ggplot object with stepwise index plots for each term added.
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
      subset(df_long, .data$term %in% prev_terms)[, c("level", "index", "term", "label")]
    }))
    expanded$shown_term <- factor(rep(term_levels, sapply(seq_along(term_levels), function(i) sum(df_long$term %in% term_levels[1:i]))),
      levels = term_levels
    )
    expanded$shown_label <- factor(label_map[as.character(expanded$shown_term)], levels = label_map[term_levels])
    ggplot(expanded, aes(x = .data$level, y = .data$index, group = .data$term, colour = .data$term)) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey") +
      geom_line() +
      geom_point() +
      facet_wrap(~shown_label, ncol = 1, scales = "fixed") +
      labs(x = obj$focus, y = "Index", colour = "Step") +
      theme(strip.text = element_text(hjust = 0)) +
      ylim(0, NA)
  } else {
    ggplot(df_long, aes(x = .data$level, y = .data$index, group = 1)) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey") +
      geom_line(color = "royalblue", alpha = 0.5) +
      geom_point(color = "royalblue") +
      facet_wrap(~label, ncol = 1, scales = "fixed") +
      labs(x = obj$focus, y = "Index") +
      theme(strip.text = element_text(hjust = 0)) +
      ylim(0, NA)
  }
}
