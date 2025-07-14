#' Check if Term is Tensor Product
#'
#' @param term_label Term label
#' @return Logical
#' @keywords internal
#' 
is_tensor_product <- function(term_label) {
  grepl("^(te|ti)\\(", term_label)
}
