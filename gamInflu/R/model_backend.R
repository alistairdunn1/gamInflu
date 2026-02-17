# model_backend.R
# Model abstraction layer for gamInflu
# Supports mgcv::gam, stats::glm, and glmmTMB model backends
#
# All functions are S3 generics with methods for "gam", "glm", and "glmmTMB" classes.
# Any file that directly accesses model internals should instead call these generics.

#' @importFrom stats sigma family model.frame coef vcov residuals fitted predict update formula logLik deviance
#' @noRd
NULL

# ============================================================================
# Backend Detection
# ============================================================================

#' @title Detect Model Backend
#' @description Determine whether a model is from mgcv, stats, or glmmTMB.
#' @param model A fitted model object
#' @return Character string: "gam" or "glmmTMB"
#' @noRd
detect_backend <- function(model) {
  if (inherits(model, "glmmTMB")) return("glmmTMB")
  if (inherits(model, "gam")) return("gam")
  if (inherits(model, "glm")) return("gam") # stats::glm treated as gam-compatible
  stop("Unsupported model class: ", paste(class(model), collapse = ", "),
       ". Supported classes: gam (mgcv), glm (stats), glmmTMB.",
       call. = FALSE)
}

#' @title Check if Model is Supported
#' @description Check whether a model object is from a supported backend.
#' @param model A fitted model object
#' @return Logical TRUE if supported
#' @noRd
is_supported_model <- function(model) {
  inherits(model, "gam") || inherits(model, "glm") || inherits(model, "glmmTMB")
}

# ============================================================================
# Model Family
# ============================================================================

#' @title Get Model Family Information
#' @description Extract family name, link function, and family object from a model.
#' @param model A fitted model object
#' @return A list with components: family (character), link (character), family_object
#' @noRd
get_model_family <- function(model) {
  UseMethod("get_model_family")
}

#' @export
get_model_family.gam <- function(model) {
  fam <- model$family
  list(
    family        = fam$family,
    link          = fam$link,
    family_object = fam
  )
}

#' @export
get_model_family.glm <- get_model_family.gam

#' @export
get_model_family.glmmTMB <- function(model) {
  fam <- family(model)
  list(
    family        = fam$family,
    link          = fam$link,
    family_object = fam
  )
}

# ============================================================================
# Model Data
# ============================================================================

#' @title Get Model Data
#' @description Extract the data frame used to fit the model.
#' @param model A fitted model object
#' @return A data.frame
#' @noRd
get_model_data <- function(model) {
  UseMethod("get_model_data")
}

#' @export
get_model_data.gam <- function(model) {
  if (!is.null(model$data)) return(as.data.frame(model$data))
  # Fallback 1: evaluate from call
  d <- tryCatch(
    eval(model$call$data, environment(formula(model))),
    error = function(e) NULL
  )
  if (!is.null(d)) return(as.data.frame(d))
  # Fallback 2: model frame (contains response + predictor columns)
  d <- tryCatch(model.frame(model), error = function(e) NULL)
  if (!is.null(d)) return(as.data.frame(d))
  # Fallback 3: reconstruct from model$model
  if (!is.null(model$model)) return(as.data.frame(model$model))
  stop("Cannot extract data from gam model.", call. = FALSE)
}

#' @export
get_model_data.glm <- get_model_data.gam

#' @export
get_model_data.glmmTMB <- function(model) {
  # glmmTMB stores data in model$frame or we can use model.frame()
  d <- tryCatch(model.frame(model), error = function(e) NULL)
  if (!is.null(d)) return(as.data.frame(d))
  # Fallback
  tryCatch(
    eval(model$call$data, environment(formula(model))),
    error = function(e) stop("Cannot extract data from glmmTMB model.", call. = FALSE)
  )
}

# ============================================================================
# Predictions
# ============================================================================

#' @title Predict on Response Scale
#' @description Get predictions on the response scale with standard errors.
#' @param model A fitted model object
#' @param newdata Data frame for predictions (uses model data if NULL)
#' @return A list with components: fit (numeric vector), se.fit (numeric vector)
#' @noRd
predict_response <- function(model, newdata = NULL) {
  UseMethod("predict_response")
}

#' @export
predict_response.gam <- function(model, newdata = NULL) {
  args <- list(model, type = "response", se.fit = TRUE)
  if (!is.null(newdata)) args$newdata <- newdata
  p <- do.call(stats::predict, args)
  list(fit = as.numeric(p$fit), se.fit = as.numeric(p$se.fit))
}

#' @export
predict_response.glm <- predict_response.gam

#' @export
predict_response.glmmTMB <- function(model, newdata = NULL) {
  args <- list(model, type = "response", se.fit = TRUE, re.form = NA)
  if (!is.null(newdata)) args$newdata <- newdata
  p <- do.call(stats::predict, args)
  se <- attr(p, "se.fit")
  if (is.null(se)) se <- rep(NA_real_, length(p))
  list(fit = as.numeric(p), se.fit = as.numeric(se))
}

#' @title Predict on Link Scale
#' @description Get predictions on the link (linear predictor) scale with standard errors.
#' @param model A fitted model object
#' @param newdata Data frame for predictions (uses model data if NULL)
#' @return A list with components: fit (numeric vector), se.fit (numeric vector)
#' @noRd
predict_link <- function(model, newdata = NULL) {
  UseMethod("predict_link")
}

#' @export
predict_link.gam <- function(model, newdata = NULL) {

  args <- list(model, type = "link", se.fit = TRUE)
  if (!is.null(newdata)) args$newdata <- newdata
  p <- do.call(stats::predict, args)
  list(fit = as.numeric(p$fit), se.fit = as.numeric(p$se.fit))
}

#' @export
predict_link.glm <- predict_link.gam

#' @export
predict_link.glmmTMB <- function(model, newdata = NULL) {
  args <- list(model, type = "link", se.fit = TRUE, re.form = NA)
  if (!is.null(newdata)) args$newdata <- newdata
  p <- do.call(stats::predict, args)
  se <- attr(p, "se.fit")
  if (is.null(se)) se <- rep(NA_real_, length(p))
  list(fit = as.numeric(p), se.fit = as.numeric(se))
}

#' @title Predict Terms (Partial Effects)
#' @description Get partial (term-by-term) predictions. This is the critical function
#' for influence analysis. mgcv supports type="terms" natively, while glmmTMB
#' must reconstruct partial predictions from the model matrix.
#' @param model A fitted model object
#' @param newdata Data frame for predictions
#' @return A list with components:
#'   fit - matrix with columns named by term, rows by observation
#'   se.fit - matrix of standard errors with same structure
#' @noRd
predict_terms <- function(model, newdata = NULL) {
  UseMethod("predict_terms")
}

#' @export
predict_terms.gam <- function(model, newdata = NULL) {
  args <- list(model, type = "terms", se.fit = TRUE)
  if (!is.null(newdata)) args$newdata <- newdata
  p <- do.call(stats::predict, args)
  list(fit = p$fit, se.fit = p$se.fit)
}

#' @export
predict_terms.glm <- predict_terms.gam

#' @export
predict_terms.glmmTMB <- function(model, newdata = NULL) {
  # glmmTMB does not support type="terms" -- reconstruct from model matrix
  # Get the fixed-effect design matrix
  frm <- formula(model)
  # Remove random effects from formula to build fixed-effects-only design matrix
  fixed_terms <- stats::delete.response(stats::terms(frm))

  if (is.null(newdata)) {
    mm <- stats::model.matrix(fixed_terms, data = model.frame(model))
  } else {
    mm <- stats::model.matrix(fixed_terms, data = newdata)
  }

  # Get fixed-effect coefficients
  beta <- glmmTMB::fixef(model)$cond

  # Map columns to term indices using the 'assign' attribute
  assign_vec <- attr(mm, "assign")
  # assign_vec: 0 = intercept, 1 = first term, 2 = second term, etc.
  term_labels <- attr(stats::terms(frm), "term.labels")

  # Number of terms (excluding intercept)
  n_terms <- length(term_labels)
  n_obs   <- nrow(mm)

  # Build matrices for fit and se
  fit_mat <- matrix(0, nrow = n_obs, ncol = n_terms)
  se_mat  <- matrix(0, nrow = n_obs, ncol = n_terms)
  colnames(fit_mat) <- term_labels
  colnames(se_mat)  <- term_labels

  # Get variance-covariance matrix for fixed effects
  V <- stats::vcov(model)$cond

  for (j in seq_len(n_terms)) {
    cols <- which(assign_vec == j)
    if (length(cols) == 0) next

    # Partial prediction for this term
    X_j <- mm[, cols, drop = FALSE]
    b_j <- beta[cols]
    partial <- as.numeric(X_j %*% b_j)

    # Centre the partial prediction (match mgcv behaviour)
    partial_centred <- partial - mean(partial)
    fit_mat[, j] <- partial_centred

    # Standard errors via sandwich: sqrt(diag(X_j %*% V_j %*% t(X_j)))
    V_j <- V[cols, cols, drop = FALSE]
    se_mat[, j] <- sqrt(pmax(0, rowSums((X_j %*% V_j) * X_j)))
  }

  list(fit = fit_mat, se.fit = se_mat)
}

# ============================================================================
# Coefficients and Variance-Covariance
# ============================================================================

#' @title Get Model Coefficients
#' @description Extract fixed-effect coefficients from the model.
#' @param model A fitted model object
#' @return A named numeric vector of coefficients
#' @noRd
get_coefficients <- function(model) {
  UseMethod("get_coefficients")
}

#' @export
get_coefficients.gam <- function(model) {
  stats::coef(model)
}

#' @export
get_coefficients.glm <- get_coefficients.gam

#' @export
get_coefficients.glmmTMB <- function(model) {
  glmmTMB::fixef(model)$cond
}

#' @title Get Variance-Covariance Matrix
#' @description Extract the variance-covariance matrix for fixed-effect coefficients.
#' @param model A fitted model object
#' @return A square matrix
#' @noRd
get_vcov <- function(model) {
  UseMethod("get_vcov")
}

#' @export
get_vcov.gam <- function(model) {
  stats::vcov(model)
}

#' @export
get_vcov.glm <- get_vcov.gam

#' @export
get_vcov.glmmTMB <- function(model) {
  stats::vcov(model)$cond
}

# ============================================================================
# Model Summary Statistics
# ============================================================================

#' @title Get Model Statistics
#' @description Extract R-squared, deviance explained, smooth table, and dispersion.
#' @param model A fitted model object
#' @return A list with components: r_sq, dev_expl, s_table, dispersion
#' @noRd
get_model_stats <- function(model) {
  UseMethod("get_model_stats")
}

#' @export
get_model_stats.gam <- function(model) {
  s <- summary(model)
  list(
    r_sq       = if (!is.null(s$r.sq)) s$r.sq else NA_real_,
    dev_expl   = if (!is.null(s$dev.expl)) s$dev.expl else {
      if (!is.null(model$null.deviance) && !is.null(model$deviance))
        (model$null.deviance - model$deviance) / model$null.deviance
      else NA_real_
    },
    s_table    = if (!is.null(s$s.table)) s$s.table else NULL,
    dispersion = if (!is.null(s$dispersion)) s$dispersion else NA_real_
  )
}

#' @export
get_model_stats.glm <- get_model_stats.gam

#' @export
get_model_stats.glmmTMB <- function(model) {
  # glmmTMB doesn't provide r.sq or dev.expl directly -- approximate
  null_ll <- tryCatch(
    as.numeric(logLik(stats::update(model, . ~ 1))),
    error = function(e) NA_real_
  )
  full_ll <- as.numeric(logLik(model))

  # McFadden's pseudo-R-squared
  r_sq <- if (!is.na(null_ll) && null_ll != 0) 1 - (full_ll / null_ll) else NA_real_

  # Deviance explained (approximate)
  if (!is.na(null_ll)) {
    null_dev <- -2 * null_ll
    full_dev <- -2 * full_ll
    dev_expl <- if (null_dev != 0) (null_dev - full_dev) / null_dev else NA_real_
  } else {
    dev_expl <- NA_real_
  }

  # Dispersion
  disp <- tryCatch(sigma(model)^2, error = function(e) NA_real_)

  list(
    r_sq       = r_sq,
    dev_expl   = dev_expl,
    s_table    = NULL, # glmmTMB has no smooth table
    dispersion = disp
  )
}

#' @title Get Model Fit Statistics
#' @description Extract null deviance, deviance, AIC, log-likelihood, residual df.
#' @param model A fitted model object
#' @return A list with components
#' @noRd
get_model_fit_stats <- function(model) {
  UseMethod("get_model_fit_stats")
}

#' @export
get_model_fit_stats.gam <- function(model) {
  list(
    null_deviance = model$null.deviance,
    deviance      = stats::deviance(model),
    aic           = stats::AIC(model),
    loglik        = as.numeric(stats::logLik(model)),
    df_residual   = model$df.residual
  )
}

#' @export
get_model_fit_stats.glm <- get_model_fit_stats.gam

#' @export
get_model_fit_stats.glmmTMB <- function(model) {
  list(
    null_deviance = NA_real_, # glmmTMB doesn't store null deviance directly
    deviance      = -2 * as.numeric(stats::logLik(model)),
    aic           = stats::AIC(model),
    loglik        = as.numeric(stats::logLik(model)),
    df_residual   = stats::df.residual(model)
  )
}

# ============================================================================
# Residuals, Fitted, and Response
# ============================================================================

#' @title Get Model Residuals
#' @description Extract residuals of a given type from the model.
#' @param model A fitted model object
#' @param type Character, residual type (e.g., "deviance", "pearson", "response", "working")
#' @return A numeric vector
#' @noRd
get_residuals <- function(model, type = "deviance") {
  UseMethod("get_residuals")
}

#' @export
get_residuals.gam <- function(model, type = "deviance") {
  stats::residuals(model, type = type)
}

#' @export
get_residuals.glm <- get_residuals.gam

#' @export
get_residuals.glmmTMB <- function(model, type = "deviance") {
  # glmmTMB supports "response" and "pearson" directly
  if (type %in% c("response", "pearson")) {
    return(stats::residuals(model, type = type))
  }
  # For "deviance", compute manually
  if (type == "deviance") {
    y  <- stats::model.response(stats::model.frame(model))
    mu <- stats::fitted(model)
    fam <- family(model)
    dev_resids <- tryCatch(
      fam$dev.resids(y, mu, 1),
      error = function(e) {
        # Fallback to pearson residuals
        warning("Could not compute deviance residuals for glmmTMB, falling back to pearson.")
        return(stats::residuals(model, type = "pearson")^2)
      }
    )
    return(sign(y - mu) * sqrt(pmax(0, dev_resids)))
  }
  # Default fallback for "working" etc.
  stats::residuals(model, type = "response")
}

#' @title Get Fitted Values
#' @description Extract fitted values from the model.
#' @param model A fitted model object
#' @return A numeric vector
#' @noRd
get_fitted <- function(model) {
  stats::fitted(model)
}

#' @title Get Model Response
#' @description Extract the response variable from the model.
#' @param model A fitted model object
#' @return A numeric vector
#' @noRd
get_model_response <- function(model) {
  UseMethod("get_model_response")
}

#' @export
get_model_response.gam <- function(model) {
  model$y
}

#' @export
get_model_response.glm <- get_model_response.gam

#' @export
get_model_response.glmmTMB <- function(model) {
  stats::model.response(stats::model.frame(model))
}

# ============================================================================
# Model Update and Fitting
# ============================================================================

#' @title Update Model
#' @description Update a model with a new formula (for stepwise building).
#' @param model A fitted model object
#' @param formula_change Formula or character specifying the change
#' @param data Data frame for refitting
#' @return An updated model of the same class
#' @noRd
update_model <- function(model, formula_change, data = NULL) {
  UseMethod("update_model")
}

#' @export
update_model.gam <- function(model, formula_change, data = NULL) {
  if (is.character(formula_change)) {
    formula_change <- stats::as.formula(formula_change)
  }
  args <- list(model, formula_change)
  if (!is.null(data)) args$data <- data
  do.call(stats::update, args)
}

#' @export
update_model.glm <- update_model.gam

#' @export
update_model.glmmTMB <- function(model, formula_change, data = NULL) {
  if (is.character(formula_change)) {
    formula_change <- stats::as.formula(formula_change)
  }
  args <- list(model, formula_change)
  if (!is.null(data)) args$data <- data
  do.call(stats::update, args)
}

#' @title Fit a New Model
#' @description Fit a new model using the same backend as an existing model.
#' Used when analyse_focus_by_group and similar functions need to fit fresh models.
#' @param model An existing fitted model (used to determine backend)
#' @param formula A formula for the new model
#' @param data Data frame for fitting
#' @param family Family specification
#' @return A fitted model object
#' @noRd
fit_model <- function(model, formula, data, family = NULL) {
  UseMethod("fit_model")
}

#' @export
fit_model.gam <- function(model, formula, data, family = NULL) {
  if (is.null(family)) family <- model$family
  mgcv::gam(formula = formula, data = data, family = family)
}

#' @export
fit_model.glm <- fit_model.gam

#' @export
fit_model.glmmTMB <- function(model, formula, data, family = NULL) {
  if (is.null(family)) family <- family(model)
  glmmTMB::glmmTMB(formula = formula, data = data, family = family)
}

# ============================================================================
# Random Effects
# ============================================================================

#' @title Extract Random Effects
#' @description Extract random effect coefficients and standard errors.
#' For mgcv, extracts from smooth objects with class "random.effect".
#' For glmmTMB, uses ranef() with condVar.
#' @param model A fitted model object
#' @return A named list; each element is a list with $coefficients and $std_errors
#' @noRd
extract_random_effects_generic <- function(model) {
  UseMethod("extract_random_effects_generic")
}

#' @export
extract_random_effects_generic.gam <- function(model) {
  re_list <- list()
  if (is.null(model$smooth)) return(re_list)
  for (sm in model$smooth) {
    if (inherits(sm, "random.effect")) {
      coefs <- model$coefficients[sm$first.para:sm$last.para]
      ses   <- sqrt(diag(model$Vp))[sm$first.para:sm$last.para]
      re_list[[sm$label]] <- list(coefficients = coefs, std_errors = ses)
    }
  }
  re_list
}

#' @export
extract_random_effects_generic.glm <- function(model) {
  # stats::glm doesn't have random effects
  list()
}

#' @export
extract_random_effects_generic.glmmTMB <- function(model) {
  re <- glmmTMB::ranef(model, condVar = TRUE)$cond
  re_list <- list()
  for (grp_name in names(re)) {
    grp_re <- re[[grp_name]]
    for (col_name in colnames(grp_re)) {
      label <- if (col_name == "(Intercept)") {
        paste0("s(", grp_name, ")")
      } else {
        paste0("s(", grp_name, ",", col_name, ")")
      }
      coefs <- grp_re[[col_name]]
      names(coefs) <- rownames(grp_re)

      # Extract posterior variances
      pv <- attr(grp_re, "postVar")
      if (!is.null(pv)) {
        col_idx <- which(colnames(grp_re) == col_name)
        ses <- sqrt(pv[col_idx, col_idx, ])
        names(ses) <- rownames(grp_re)
      } else {
        ses <- rep(NA_real_, length(coefs))
        names(ses) <- rownames(grp_re)
      }
      re_list[[label]] <- list(coefficients = coefs, std_errors = ses)
    }
  }
  re_list
}

# ============================================================================
# Smooth EDF
# ============================================================================

#' @title Get Smooth EDF
#' @description Get the total effective degrees of freedom from smooth terms.
#' Returns 0 for glmmTMB (which uses fixed-df spline bases).
#' @param model A fitted model object
#' @return Numeric value
#' @noRd
get_smooth_edf <- function(model) {
  UseMethod("get_smooth_edf")
}

#' @export
get_smooth_edf.gam <- function(model) {
  s_tab <- tryCatch(summary(model)$s.table, silent = TRUE)
  if (inherits(s_tab, "try-error") || is.null(s_tab)) return(0)
  if (is.matrix(s_tab) || is.data.frame(s_tab)) {
    edf_col <- if ("edf" %in% colnames(s_tab)) s_tab[, "edf"] else if (ncol(s_tab) >= 1) s_tab[, 1] else 0
    return(sum(as.numeric(edf_col), na.rm = TRUE))
  }
  0
}

#' @export
get_smooth_edf.glm <- function(model) { 0 }

#' @export
get_smooth_edf.glmmTMB <- function(model) { 0 }

# ============================================================================
# Dispersion Parameters
# ============================================================================

#' @title Get Dispersion Parameters
#' @description Extract dispersion and theta parameters from the model.
#' @param model A fitted model object
#' @return A list with components: dispersion, theta (NULL if not applicable)
#' @noRd
get_dispersion_params <- function(model) {
  UseMethod("get_dispersion_params")
}

#' @export
get_dispersion_params.gam <- function(model) {
  s <- summary(model)
  theta <- tryCatch(model$family$getTheta(TRUE), error = function(e) NULL)
  list(
    dispersion = if (!is.null(s$dispersion)) s$dispersion else NA_real_,
    theta      = theta
  )
}

#' @export
get_dispersion_params.glm <- get_dispersion_params.gam

#' @export
get_dispersion_params.glmmTMB <- function(model) {
  disp <- tryCatch(sigma(model)^2, error = function(e) NA_real_)
  # For NB models, theta is the overdispersion parameter
  theta <- tryCatch({
    fam <- family(model)
    if (grepl("nbinom", fam$family, ignore.case = TRUE)) {
      glmmTMB::getME(model, "theta")
    } else {
      NULL
    }
  }, error = function(e) NULL)

  list(dispersion = disp, theta = theta)
}

# ============================================================================
# Model Offset
# ============================================================================

#' @title Get Model Offset
#' @description Extract offset from a model.
#' @param model A fitted model object
#' @return Numeric vector or NULL
#' @noRd
get_model_offset <- function(model) {
  UseMethod("get_model_offset")
}

#' @export
get_model_offset.gam <- function(model) {
  model$offset
}

#' @export
get_model_offset.glm <- get_model_offset.gam

#' @export
get_model_offset.glmmTMB <- function(model) {
  tryCatch({
    of <- stats::model.offset(stats::model.frame(model))
    if (length(of) == 0) NULL else of
  }, error = function(e) NULL)
}
