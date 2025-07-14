#' Print Method for gam_influence Objects
#'
#' @param step_effects A list of step effects
#' @param x A gam_influence object
#' @param show_previous Logical indicating whether to show previous steps
#' @keywords internal
#'
prepare_step_data <- function(step_effects, x, show_previous = TRUE) {
  if (length(step_effects) == 0) {
    return(data.frame())
  }

  # Get clean level names
  clean_levels <- extract_clean_levels(names(x$focus_effects$relative_effects), x$focus)

  # Prepare data for all steps
  all_data <- data.frame()

  for (i in seq_along(step_effects)) {
    step <- step_effects[[i]]

    if (!is.null(step$effects) && length(step$effects) > 0) {
      # Current step data
      current_data <- data.frame(
        step_number = i,
        step_name = paste(i, ". ", step$term_added, sep = ""),
        term_added = step$term_added,
        level = clean_levels,
        level_number = seq_along(clean_levels),
        effect = step$effects,
        is_current = TRUE,
        is_previous = FALSE,
        line_type = "Current"
      )

      all_data <- rbind(all_data, current_data)

      # Add previous steps data if requested
      if (show_previous && i > 1) {
        for (j in 1:(i - 1)) {
          prev_step <- step_effects[[j]]
          if (!is.null(prev_step$effects) && length(prev_step$effects) > 0) {
            previous_data <- data.frame(
              step_number = i,
              step_name = paste(i, ". ", step$term_added, sep = ""),
              term_added = step$term_added,
              level = clean_levels,
              level_number = seq_along(clean_levels),
              effect = prev_step$effects,
              is_current = FALSE,
              is_previous = TRUE,
              line_type = paste("Step", j)
            )
            all_data <- rbind(all_data, previous_data)
          }
        }
      }
    }
  }

  return(all_data)
}
