#' @title Subplot for Focus Effects
#' @description Internal function to plot predicted effects for focus terms.
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param preds_df Data frame of predictions.
#' @param se_df Data frame of standard errors.
#' @param term_vars Variables in the term.
#' @param se_col Standard error column(s).
#' @param islog Logical; if TRUE, exponentiate effects.
#' @noRd
subplot_focus_effect <- function(obj, t, term_vars, cdi = FALSE) {
  message("Plotting focus effect for term: ", t)
  if (t != obj$focus) {
    message(paste(t, " is not a focus term", sep = ""))
    return(NULL)
  }

  islog <- isTRUE(obj$islog)
  obj_data <- obj$data
  by_var <- sub(".*by\\s*=\\s*([^,\\)]+).*", "\\1", t)
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se

  effect <- preds_df[[t]]
  se <- se_df[[t]]
  if (islog) {
    lower <- exp(effect - 1.96 * se)
    upper <- exp(effect + 1.96 * se)
    effect <- exp(effect)
    ylim <- c(0, NA)
  } else {
    lower <- effect - 1.96 * se
    upper <- effect + 1.96 * se
    ylim <- c(NA_real_, NA_real_)
  }

  # Handle potential dimension mismatch between data and predictions (subset analysis)
  if (length(effect) != nrow(obj$data)) {
    # This indicates subset analysis was performed
    message("Detected subset analysis: adjusting focus variable levels to match predictions")

    # Get all unique levels from the current data
    all_levels <- unique(obj$data[[term_vars[1]]][!is.na(obj$data[[term_vars[1]]])])

    # If we have the right number of levels, use them
    if (length(all_levels) == length(effect)) {
      df <- data.frame(level = all_levels, effect = effect, lower = lower, upper = upper)
    } else {
      # Create levels based on sorted unique values
      sorted_levels <- sort(all_levels)
      if (length(sorted_levels) >= length(effect)) {
        used_levels <- sorted_levels[seq_len(length(effect))]
        df <- data.frame(level = used_levels, effect = effect, lower = lower, upper = upper)
      } else {
        # Fallback: use the focus variable sequence if it's numeric
        if (is.numeric(all_levels)) {
          df <- data.frame(level = seq_len(length(effect)), effect = effect, lower = lower, upper = upper)
        } else {
          df <- data.frame(level = paste("Level", seq_len(length(effect))), effect = effect, lower = lower, upper = upper)
        }
        warning("Could not match focus variable levels to predictions. Using generic labels.")
      }
    }
  } else {
    # Normal case - dimensions match
    df <- data.frame(level = obj$data[[term_vars[1]]], effect = effect, lower = lower, upper = upper)
  }

  # Convert df$level to numeric if it is numeric-like (e.g., factor or character with numeric values)
  if (!is.numeric(df$level) && all(!is.na(as.numeric(as.character(df$level))))) {
    df$level <- as.numeric(as.character(df$level))
  }

  p_coef <- ggplot2::ggplot(df, ggplot2::aes(x = level, y = effect)) +
    ggplot2::geom_point(colour = "royalblue") +
    ggplot2::geom_line(colour = "royalblue") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "royalblue", alpha = 0.2) +
    ggplot2::ylim(ylim) +
    ggplot2::labs(y = "Index")
  if (cdi) {
    p_coef <- p_coef + ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    )
  } else {
    p_coef <- p_coef +
      ggplot2::xlab(term_vars[1])
  }
  return(p_coef)
}
