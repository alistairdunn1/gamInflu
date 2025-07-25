#' @title Standardisation Plot
#' @description Creates a standardisation plot comparing the unstandardised (raw) index to the final standardised index for the focus term.
#' @param obj A `gam_influence` object containing calculated indices.
#' @return A ggplot object showing both unstandardised and standardised indices, with confidence ribbon.
#' @importFrom ggplot2 ggplot aes geom_hline geom_line geom_point geom_ribbon labs scale_colour_manual scale_y_continuous
#' @importFrom rlang .data
#' @export
plot_standardisation <- function(obj) {
  df <- obj$calculated$indices
  if (is.null(df)) {
    stop("No indices calculated. Please run `calculate_influence()` first.", call. = FALSE)
  }
  if (!("unstan" %in% names(df)) || !("stan" %in% names(df))) {
    stop("Data frame must contain 'unstan' and 'stan' columns.", call. = FALSE)
  }
  if (!("level" %in% names(df))) {
    stop("Data frame must contain a 'level' column for the focus term.", call. = FALSE)
  }
  if (!("stanLower" %in% names(df)) || !("stanUpper" %in% names(df))) {
    stop("Data frame must contain 'stanLower' and 'stanUpper' columns for the standardised index range.", call. = FALSE)
  }
  if (is.null(obj$focus)) {
    stop("Focus term is not set. Please set the focus term in the gam_influence object.", call. = FALSE)
  }

  # Convert level to numeric if possible
  if (is.factor(df$level) && all(!is.na(as.numeric(as.character(levels(df$level)))))) {
    df$level <- as.numeric(as.character(df$level))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$level, group = 1)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 1), linetype = "dashed", colour = "grey") +
    ggplot2::geom_line(ggplot2::aes(y = .data$unstan, colour = "Unstandardised")) +
    ggplot2::geom_point(ggplot2::aes(y = .data$unstan, colour = "Unstandardised")) +
    ggplot2::geom_line(ggplot2::aes(y = .data$stan, colour = "Standardised")) +
    ggplot2::geom_point(ggplot2::aes(y = .data$stan, colour = "Standardised")) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$stanLower, ymax = .data$stanUpper), fill = "royalblue", alpha = 0.2) +
    ggplot2::labs(x = obj$focus, y = "Index", colour = "") +
    ggplot2::scale_colour_manual(values = c("Unstandardised" = "grey40", "Standardised" = "royalblue")) +
    ggplot2::scale_y_continuous(limits = c(0, NA))

  return(p)
}
