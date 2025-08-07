#' @title QQ Plot for Random Effects
#' @description QQ plot for random effect coefficients in a gam_influence object.
#' @param obj A `gam_influence` object.
#' @param term Random effect term name.
#' @return A ggplot object showing the QQ plot of random effect coefficients.
#' @importFrom ggplot2 aes geom_point geom_abline labs geom_histogram geom_density geom_vline geom_errorbar geom_hline ylim after_stat
#' @importFrom patchwork plot_spacer
#' @importFrom stats qnorm ppoints sd
#' @noRd
subplot_random_effect_qq <- function(obj, term) {
  message("Plotting QQ plot for random effect term: ", term)
  re_list <- extract_random_effects(obj)
  re_label <- .match_re_label(term, re_list)
  if (is.null(re_label)) {
    message("No random effect found for term: ", term)
    return(patchwork::plot_spacer())
  }
  re_info <- re_list[[re_label]]
  coeffs <- re_info$coefficients
  coeffs_std <- as.numeric(scale(coeffs))
  n <- length(coeffs_std)
  theoretical <- qnorm(ppoints(n))
  df <- data.frame(theoretical = sort(theoretical), sample = sort(coeffs_std))
  p_re <- ggplot(df, aes(x = .data$theoretical, y = .data$sample)) +
    geom_point(colour = "royalblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    labs(x = "Theoretical Quantiles", y = "Sample quantiles")
  return(p_re)
}

#' @title Histogram for Random Effects
#' @description Histogram/density plot for random effect coefficients in a gam_influence object.
#' @param obj A `gam_influence` object.
#' @param term Random effect term name.
#' @noRd
subplot_random_effect_histogram <- function(obj, term) {
  message("Plotting histogram for random effect term: ", term)
  re_list <- extract_random_effects(obj)
  re_label <- .match_re_label(term, re_list)
  if (is.null(re_label)) {
    message("No random effect found for term: ", term)
    return(patchwork::plot_spacer())
  }
  re_info <- re_list[[re_label]]
  coeffs <- re_info$coefficients
  df <- data.frame(coefficient = as.numeric(coeffs))
  p_re <- ggplot(df, aes(x = .data$coefficient)) +
    geom_histogram(aes(y = after_stat(density)), bins = 15, alpha = 0.7, fill = "royalblue") +
    geom_density(colour = "royalblue") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(x = "Random Effect Coefficient", y = "Density")
  return(p_re)
}

#' @title Points for Random Effects
#' @description Points plot for random effect coefficients in a gam_influence object.
#' @param obj A `gam_influence` object.
#' @param term Random effect term name.
#' @noRd
subplot_random_effect_points <- function(obj, term) {
  message("Plotting points for random effect term: ", term)
  re_list <- extract_random_effects(obj)
  re_label <- .match_re_label(term, re_list)
  if (is.null(re_label)) {
    message("No random effect found for term: ", term)
    return(patchwork::plot_spacer())
  }
  re_info <- re_list[[re_label]]
  coeffs <- re_info$coefficients
  df <- data.frame(coefficient = as.numeric(coeffs), ID = seq_along(coeffs))

  if (obj$islog) {
    df$lower <- exp(df$coefficient - 1.96 * re_info$std_errors)
    df$upper <- exp(df$coefficient + 1.96 * re_info$std_errors)
    df$coefficient <- exp(df$coefficient)
    ylim <- c(0, NA)
  } else {
    df$lower <- df$coefficient - 1.96 * re_info$std_errors
    df$upper <- df$coefficient + 1.96 * re_info$std_errors
    ylim <- c(NA_real_, NA_real_)
  }

  p_re <- ggplot(df, aes(x = ID, y = coefficient)) +
    geom_point(colour = "royalblue") +
    geom_errorbar(aes(ymin = lower, ymax = upper), colour = "royalblue", alpha = 0.5, width = 0.2) +
    labs(x = "Level", y = "Random effect") +
    ylim(ylim)

  return(p_re)
}

#' @title Caterpillar Plot for Random Effects
#' @description Caterpillar plot for random effect coefficients in a gam_influence object.
#' @param obj A `gam_influence` object.
#' @param term Random effect term name.
#' @param conf_level Confidence level for error bars.
#' @noRd
subplot_random_effect_caterpillar <- function(obj, term, conf_level = 0.95) {
  message("Plotting caterpillar plot for random effect term: ", term)
  re_list <- extract_random_effects(obj)
  re_label <- .match_re_label(term, re_list)
  if (is.null(re_label)) {
    message("No random effect found for term: ", term)
    return(patchwork::plot_spacer())
  }
  re_info <- re_list[[re_label]]
  coeffs <- re_info$coefficients
  ses <- re_info$std_errors
  alpha <- 1 - conf_level
  z_score <- qnorm(1 - alpha / 2)
  df <- data.frame(
    level = names(coeffs),
    estimate = as.numeric(coeffs),
    lower = as.numeric(coeffs) - z_score * as.numeric(ses),
    upper = as.numeric(coeffs) + z_score * as.numeric(ses)
  )
  df$level <- factor(df$level, levels = df$level)
  p_re <- ggplot(df, aes(x = level, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = lower, ymax = upper), colour = "royalblue", alpha = 0.5, width = 0.3) +
    geom_point(colour = "royalblue") +
    labs(x = "Level", y = "Random effect")
  return(p_re)
}

#' @title Extract Random Effects from mgcv Model
#' @description Extracts random effect coefficients and standard errors from a fitted mgcv model.
#' @param obj A `gam_influence` object.
#' @return A named list of random effect coefficients and standard errors.
#' @noRd
extract_random_effects <- function(obj) {
  model <- obj$model
  re_list <- list()
  for (sm in model$smooth) {
    if (inherits(sm, "random.effect")) {
      coefs <- model$coefficients[sm$first.para:sm$last.para]
      ses <- sqrt(diag(model$Vp))[sm$first.para:sm$last.para]
      re_list[[sm$label]] <- list(coefficients = coefs, std_errors = ses)
    }
  }
  return(re_list)
}

#' @title Match Random Effect Label
#' @description Helper function to match a term name to a random effect label in the extracted random effects list.
#' @param term Character string of the term name to match.
#' @param re_list Named list of random effects from extract_random_effects().
#' @return Character string of the matched label, or NULL if no match found.
#' @noRd
.match_re_label <- function(term, re_list) {
  idx <- which(vapply(names(re_list), function(nm) grepl(term, nm, fixed = TRUE), logical(1)))
  if (length(idx) == 0) {
    return(NULL)
  }
  names(re_list)[idx[1]]
}
