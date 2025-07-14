#' Check if term is related to the focus term
#'
#' @param term_name Name of the term to check
#' @param focus_term Name of the focus term
#' @return Logical indicating if term is focus-related
#' @keywords internal
#' 
is_focus_related <- function(term_name, focus_term) {
  # Check if the term name contains the focus term
  grepl(focus_term, term_name, fixed = TRUE)
}
