#' @title Random Effects Diagnostics
#' @description Plot diagnostics for random effects in a gam_influence object.
#' @param obj A `gam_influence` object.
#' @param term Random effect term name.
#' @param re_type Type of diagnostic: "points", "qq", "hist", or "caterpillar"
#' @export
plot_re <- function(obj, term, re_type = "qq") {
  # --- Setup ---
  if (is.numeric(term) && length(term) == 1 && term == as.integer(term)) {
    all_terms <- get_terms(obj, full = TRUE)
    if (term > length(all_terms) || term < 1) {
      stop("Term index out of bounds. There are ", length(all_terms), " valid terms:\n  ", paste(all_terms, collapse = "\n  "), call. = FALSE)
    }
    term <- all_terms[term]
  }

  re_type <- pmatch(re_type, c("points", "qq", "hist", "caterpillar"), nomatch = NA)
  if (is.na(re_type)) {
    message("Unsupported re_type for random effect plot. Use 'points', 'qq', 'hist', or 'caterpillar'")
    p_coef <- plot_spacer()
  } else if (re_type == 1) {
    p_coef <- subplot_random_effect_points(obj, term_vars)
  } else if (re_type == 2) {
    p_coef <- subplot_random_effect_qq(obj, term_vars)
  } else if (re_type == 3) {
    p_coef <- subplot_random_effect_histogram(obj, term_vars)
  } else if (re_type == 4) {
    p_coef <- subplot_random_effect_caterpillar(obj, term_vars)
  }
}
