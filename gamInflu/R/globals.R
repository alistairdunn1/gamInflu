# Globals for R CMD check compliance
# Variables used in ggplot2 aesthetics and dplyr operations
#' @noRd
utils::globalVariables(c(
  # Variable names used in ggplot2 aesthetics
  "level", "index", "term", "influence", "effect", "x", "y", "z",
  "lower", "upper", "var1", "var2", "estimate", "density", "ID", "coefficient",

  # Variables used in data manipulation
  "term_vars", "xx"
))
