#' @title Standardisation Plot
#' @description Creates a standardisation plot comparing the unstandardised (raw) index to the final standardised index for the focus term.
#' @param obj A `gam_influence` object containing calculated indices.
#' @return A ggplot object showing both unstandardised and standardised indices, with confidence ribbon.
#' @export
#' @describeIn plot.gam_influence Creates a standardisation plot.
plot_standardisation <- function(obj) {
  df <- obj$calculated$indices
  if (is.null(df)) {
    stop("No indices calculated. Please run `calculate_influence()` first.")
  }
  if (!("unstan" %in% names(df)) || !("stan" %in% names(df))) {
    stop("Data frame must contain 'unstan' and 'stan' columns.")
  }
  if (!("level" %in% names(df))) {
    stop("Data frame must contain a 'level' column for the focus term.")
  }
  if (!("stanLower" %in% names(df)) || !("stanUpper" %in% names(df))) {
    stop("Data frame must contain 'stanLower' and 'stanUpper' columns for the standardised index range.")
  }
  if (is.null(obj$focus)) {
    stop("Focus term is not set. Please set the focus term in the gam_influence object.")
  }

  # Convert level to numeric if possible
  if (is.factor(df$level) && all(!is.na(as.numeric(as.character(levels(df$level)))))) {
    df$level <- as.numeric(as.character(df$level))
  }

  ggplot(df, aes(x = level, group = 1)) +
    geom_hline(aes(yintercept = 1), linetype = "dashed", colour = "grey") +
    geom_line(aes(y = unstan, colour = "Unstandardised")) +
    geom_point(aes(y = unstan, colour = "Unstandardised")) +
    geom_line(aes(y = stan, colour = "Standardised")) +
    geom_point(aes(y = stan, colour = "Standardised")) +
    geom_ribbon(aes(ymin = stanLower, ymax = stanUpper), alpha = 0.2, fill = "royalblue") +
    labs(x = obj$focus, y = "Index", colour = "") +
    scale_y_continuous(limits = c(0, NA))
}
