#' @title Subplot for Influence Plot
#' @description Internal function to plot the influence for a term in CDI plots.
#' @param influ_data Data frame with influence values.
#' @return A ggplot object.
#' @noRd
subplot_influence <- function(obj, term, focus_var, cdi = FALSE) {
  message("Generating influence plot for term: ", term)

  # --- Data Preparation ---
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se
  obj_data <- obj$data
  focus_var <- obj$focus
  model_terms <- get_terms(obj, full = FALSE)

  # Identify variables associated with the term (can be 1 or 2 for smooths)
  term_vars <- all.vars(rlang::parse_expr(term))
  term_vars_in_model <- term_vars[term_vars %in% model_terms]
  if (length(term_vars_in_model) == 0) {
    stop(paste0(
      "None of the variables in the supplied term ('", term, "') match the model term.\n",
      "Valid terms are: ", paste(setdiff(model_terms, obj$focus), collapse = ", ")
    ))
  }
  if (any(term_vars_in_model == focus_var)) {
    stop("CDI plot cannot be generated for the focus term itself.")
  }

  influ_data <- subset(obj$calculated$influences, term %in% term_vars_in_model) %>%
    dplyr::mutate(influence = exp(influence)) %>%
    dplyr::mutate(influence = influence / mean(influence))

  xx <<- influ_data


  ylim <- c(
    pmin(0.75, min(influ_data$influence, na.rm = TRUE)),
    pmax(1.25, max(influ_data$influence, na.rm = TRUE))
  )

  if (length(unique(influ_data$term)) == 1) {
    ggplot(influ_data, aes(x = level, y = influence, group = term)) +
      geom_line(colour = "royalblue", alpha = 0.5) +
      geom_point(colour = "royalblue") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Influence") +
      theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
      ) +
      ylim(ylim) +
      coord_flip()
  } else {
    ggplot(influ_data, aes(x = level, y = influence, group = term, colour = term)) +
      geom_line(alpha = 0.5) +
      geom_point() +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Influence", colour = "Term") +
      theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
      ) +
      ylim(ylim) +
      coord_flip() +
      theme(
        legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = rel(0.8))
      )
  }
}
