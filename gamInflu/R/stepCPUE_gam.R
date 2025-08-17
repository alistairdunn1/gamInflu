#' Stepwise Model Selection for mgcv GAM Objects Using R-squared Criterion
#'
#' @description
#' Performs stepwise model selection for mgcv GAM objects using R-squared change as the criterion.
#' This function is adapted for GAM objects and uses R-squared change rather than AIC for model selection.
#' It can perform forward, backward, or bidirectional stepwise selection.
#'
#' @param object An mgcv GAM object (fitted with \code{\link[mgcv]{gam}}).
#' @param scope Either a list with components \code{upper} and \code{lower} specifying the model scope,
#'   or a single formula specifying the upper scope. If missing, backward selection only is performed.
#' @param r2.change Numeric. The minimum change in R-squared required to add or drop a term.
#'   Default is 0.005.
#' @param scale Not used in GAM context, kept for compatibility. Default is 0.
#' @param direction Character string specifying the direction of selection:
#'   "both" (default), "backward", or "forward".
#' @param trace Numeric. Controls the amount of output during selection:
#'   \itemize{
#'     \item 0: No output
#'     \item 1: Basic progress output (default)
#'     \item 2: Detailed progress with R-squared values
#'     \item 3: Very detailed output with individual model fits
#'   }
#' @param keep A function specifying what information to keep from each fitted model.
#'   Default is \code{deviance}.
#' @param steps Maximum number of selection steps to perform. Default is 1000.
#' @param ... Additional arguments passed to \code{\link[mgcv]{update.gam}}.
#'
#' @details
#' This function adapts the stepwise selection procedure for mgcv GAM objects.
#' Unlike traditional AIC-based selection, this uses R-squared change as the criterion,
#' which can be more appropriate for GAMs where the effective degrees of freedom
#' may be difficult to determine precisely.
#'
#' The R-squared is calculated as: (null.deviance - deviance) / null.deviance
#'
#' For GAM objects, the function handles:
#' \itemize{
#'   \item Smooth terms (s(), te(), ti(), etc.)
#'   \item Parametric terms
#'   \item Random effects (s(..., bs="re"))
#'   \item Tensor products and interactions
#' }
#'
#' @return
#' An mgcv GAM object representing the final selected model, with an additional
#' \code{anova} component containing the stepwise selection path showing:
#' \itemize{
#'   \item Step: The term added or removed
#'   \item Df: Change in degrees of freedom
#'   \item Deviance: Change in deviance
#'   \item Resid. Df: Residual degrees of freedom
#'   \item Resid. Dev: Residual deviance
#'   \item r.squared: Model R-squared at each step
#'   \item aic: AIC value at each step
#' }
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#'
#' # Create example data
#' set.seed(123)
#' n <- 200
#' dat <- data.frame(
#'   x1 = runif(n, 0, 10),
#'   x2 = runif(n, 0, 10),
#'   x3 = runif(n, 0, 10),
#'   f1 = factor(sample(1:3, n, replace = TRUE))
#' )
#' dat$y <- 2 + sin(dat$x1) + 0.5 * dat$x2 + as.numeric(dat$f1) + rnorm(n, 0, 0.5)
#'
#' # Fit initial model
#' m1 <- gam(y ~ s(x1), data = dat)
#'
#' # Define scope for stepwise selection
#' scope <- list(
#'   lower = ~ s(x1),
#'   upper = ~ s(x1) + s(x2) + s(x3) + f1
#' )
#'
#' # Perform stepwise selection
#' m_step <- stepCPUE_gam(m1, scope = scope, r2.change = 0.01, trace = 2)
#'
#' # View selection path
#' print(m_step$anova)
#' }
#'
#' @seealso
#' \code{\link[mgcv]{gam}}, \code{\link[mgcv]{update.gam}}, \code{\link[stats]{step}}
#'
#' @importFrom mgcv gam
#' @importFrom stats AIC deviance formula terms update
#' @export
stepCPUE_gam <- function(object, scope, r2.change = 0.005, scale = 0,
                         direction = c("both", "backward", "forward"),
                         trace = 1, keep = deviance, steps = 1000, ...) {
  # Check if object is a GAM
  if (!inherits(object, "gam")) {
    stop("Object must be an mgcv GAM object fitted with gam()")
  }

  # Internal function to extract R-squared and other statistics from GAM
  extractCPUE_gam <- function(fit, scale = 0, ...) {
    n <- length(fit$residuals)
    # For GAMs, use the summary to get effective degrees of freedom
    edf <- sum(fit$edf) + length(fit$coefficients) - sum(fit$edf)

    # Calculate R-squared for GAMs
    if (!is.null(fit$null.deviance) && !is.null(fit$deviance)) {
      r.squared <- 1 - (fit$deviance / fit$null.deviance)
    } else {
      # Alternative calculation for GAMs
      r.squared <- summary(fit)$r.sq
    }

    c(edf, r.squared)
  }

  # Internal function to calculate R-squared for GAMs
  r.squared_gam <- function(GAM, decimals = 3) {
    if (!is.null(GAM$null.deviance) && !is.null(GAM$deviance)) {
      RES <- round(1 - (GAM$deviance / GAM$null.deviance), decimals)
    } else {
      # Use summary r.sq for GAMs
      RES <- round(summary(GAM)$r.sq, decimals)
    }
    return(RES)
  }

  # Function for backward selection (dropping terms)
  AIC.drop1_gam <- function(fit, Terms, scale, trace, ...) {
    n <- length(Terms)
    ans <- matrix(nrow = n + 1, ncol = 2)
    dimnames(ans) <- list(c("<none>", paste("-", Terms, sep = "")), c("df", "R.sq"))
    ans[1, ] <- extractCPUE_gam(fit, scale, ...)

    if (n == 0) {
      return(ans)
    }

    i <- 2
    for (tt in Terms) {
      if (trace > 1) {
        cat("trying -", tt, "\n")
      } else if (trace > 0) {
        cat(".")
      }

      # Try to update the GAM by removing the term
      tryCatch(
        {
          nfit <- update(fit, paste("~ . -", tt), ...)
          ans[i, ] <- extractCPUE_gam(nfit, scale, ...)
        },
        error = function(e) {
          # If update fails, set to NA
          ans[i, ] <- c(NA, NA)
          if (trace > 0) {
            cat("Failed to remove", tt, "\n")
          }
        }
      )

      if (trace > 2) {
        print(ans[i, ])
      }
      i <- i + 1
    }
    ans
  }

  # Function for forward selection (adding terms)
  AIC.add1_gam <- function(fit, Terms, scale, trace, ...) {
    n <- length(Terms)
    ans <- matrix(nrow = n + 1, ncol = 2)
    t2 <- if (length(Terms)) paste("+", Terms, sep = "") else NULL
    dimnames(ans) <- list(c("<none>", t2), c("df", "R.sq"))
    ans[1, ] <- extractCPUE_gam(fit, scale, ...)

    if (n == 0) {
      return(ans)
    }

    i <- 2
    for (tt in Terms) {
      if (trace > 1) {
        cat("trying +", tt, "\n")
      } else if (trace > 0) {
        cat(".")
      }

      # Try to update the GAM by adding the term
      tryCatch(
        {
          nfit <- update(fit, paste("~ . +", tt), ...)
          ans[i, ] <- extractCPUE_gam(nfit, scale, ...)
        },
        error = function(e) {
          # If update fails, set to NA
          ans[i, ] <- c(NA, NA)
          if (trace > 0) {
            cat("Failed to add", tt, "\n")
          }
        }
      )

      if (trace > 2) {
        print(ans[i, ])
      }
      i <- i + 1
    }
    ans
  }

  # Function to rearrange keep information
  re.arrange <- function(keep) {
    if (is.null(keep) || length(keep) == 0) {
      return(NULL)
    }
    namr <- names(k1 <- keep[[1]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, namc))
  }

  # Function to create the step summary
  make.step_gam <- function(models, fit, object) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, diff(rd))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, diff(rdf))
    r2_vals <- sapply(models, "[[", "R.sq")
    aic_vals <- sapply(models, "[[", "aic")

    if (trace > 0) {
      cat("\nModel selection path:\n")
      print(data.frame(R.squared = r2_vals, AIC = aic_vals))
    }

    heading <- c(
      "Stepwise Model Path \nAnalysis of Deviance Table",
      "\nInitial Model:", deparse(as.vector(formula(object))),
      "\nFinal Model:", deparse(as.vector(formula(fit))), "\n"
    )

    aod <- data.frame(
      Step = change,
      Df = ddf,
      Deviance = dd,
      "Resid. Df" = rdf,
      "Resid. Dev" = rd,
      r.squared = r2_vals,
      aic = aic_vals,
      check.names = FALSE
    )

    attr(aod, "heading") <- heading
    attr(aod, "class") <- c("anova", "data.frame")
    fit$anova <- aod
    fit
  }

  # Parse direction argument
  if (missing(direction)) {
    direction <- "both"
  } else {
    direction <- match.arg(direction)
  }

  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"

  # Parse scope argument
  if (missing(scope)) {
    fdrop <- numeric(0)
    fadd <- NULL
  } else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower)) {
        attr(terms(update.formula(object, fdrop)), "term.labels")
      } else {
        character(0)
      }
      fadd <- if (!is.null(fadd <- scope$upper)) {
        attr(terms(update.formula(object, fadd)), "term.labels")
      } else {
        NULL
      }
    } else {
      fadd <- if (!is.null(fadd <- scope)) {
        attr(terms(update.formula(object, scope)), "term.labels")
      } else {
        NULL
      }
      fdrop <- character(0)
    }
  }

  if (is.null(fadd)) {
    backward <- TRUE
    forward <- FALSE
  }

  # Initialize tracking variables
  models <- vector("list", steps)
  if (!is.null(keep)) {
    keep.list <- vector("list", steps)
  }

  n <- length(object$residuals)
  fit <- object

  # Get initial statistics
  initial_stats <- extractCPUE_gam(fit, scale, ...)
  edf <- initial_stats[1]
  bR2 <- initial_stats[2]
  aic <- AIC(fit)

  nm <- 1
  Terms <- fit$terms

  if (trace > 0) {
    cat("\n\nTesting for a change in R-squared of less than", format(round(r2.change, 4)), "\n")
  }

  # Store initial model
  models[[nm]] <- list(
    deviance = deviance(fit),
    df.resid = n - edf,
    change = "",
    R.sq = bR2,
    aic = aic
  )

  if (!is.null(keep)) {
    keep.list[[nm]] <- keep(fit, bR2)
  }

  count.steps <- 0
  current_terms <- attr(Terms, "term.labels")

  # Main stepwise loop
  while (steps > 0) {
    steps <- steps - 1
    R2_old <- bR2
    bfit <- fit

    # Get current model scope
    current_terms <- attr(fit$terms, "term.labels")

    # Determine which terms can be dropped (not in fdrop)
    droppable <- setdiff(current_terms, fdrop)

    # Determine which terms can be added (in fadd but not current)
    if (!is.null(fadd)) {
      addable <- setdiff(fadd, current_terms)
    } else {
      addable <- character(0)
    }

    aod.drop <- NULL
    aod.add <- NULL
    aod <- NULL
    change <- NULL

    # Try backward step
    if (backward && length(droppable) > 0) {
      aod.drop <- AIC.drop1_gam(fit, droppable, scale = scale, trace = trace, ...)
    }

    # Try forward step
    if (forward && length(addable) > 0) {
      aod.add <- AIC.add1_gam(fit, addable, scale = scale, trace = trace, ...)
    }

    # If no moves possible, break
    if (is.null(aod.drop) && is.null(aod.add)) {
      break
    }

    # Check for beneficial backward step
    if (!is.null(aod.drop) && nrow(aod.drop) > 1) {
      # Find best improvement (highest R-squared after removal)
      r2_gains <- aod.drop[-1, "R.sq"] - aod.drop[1, "R.sq"]
      if (any(!is.na(r2_gains)) && max(r2_gains, na.rm = TRUE) >= -r2.change) {
        best_drop <- which.max(r2_gains) + 1
        change <- dimnames(aod.drop)[[1]][best_drop]
        if (trace > 0) {
          cat(paste("\nDrop term", change, "\n"))
        }
        aod <- aod.drop
      }
    }

    # Check for beneficial forward step (only if no drop was selected)
    if (is.null(change) && !is.null(aod.add) && nrow(aod.add) > 1) {
      # Find best improvement (highest R-squared after addition)
      r2_gains <- aod.add[-1, "R.sq"] - aod.add[1, "R.sq"]
      if (any(!is.na(r2_gains)) && max(r2_gains, na.rm = TRUE) >= r2.change) {
        best_add <- which.max(r2_gains) + 1
        change <- dimnames(aod.add)[[1]][best_add]
        if (trace > 0) {
          cat(paste("\nAdd term", change, "\n"))
        }
        aod <- aod.add
      }
    }

    # If no beneficial change found, stop
    if (is.null(change)) {
      if (trace > 0) {
        cat("\nNo beneficial changes found. Stopping.\n")
        if (!is.null(aod.drop)) {
          cat("Drop options:\n")
          print(data.frame(df = aod.drop[, 1], R.squared = aod.drop[, 2]))
        }
        if (!is.null(aod.add)) {
          cat("Add options:\n")
          print(data.frame(df = aod.add[, 1], R.squared = aod.add[, 2]))
        }
      }
      break
    }

    # Apply the change
    if (trace > 1) {
      if (!is.null(aod)) {
        print(data.frame(df = aod[, 1], R.squared = aod[, 2]))
      }
    }

    # Update the model
    tryCatch(
      {
        fit <- update(fit, eval(parse(text = paste("~ .", change))), ...)
        Terms <- fit$terms

        # Get new statistics
        new_stats <- extractCPUE_gam(fit, scale, ...)
        edf <- new_stats[1]
        bR2 <- new_stats[2]
        aic <- AIC(fit)

        # Store the new model
        nm <- nm + 1
        models[[nm]] <- list(
          deviance = deviance(fit),
          df.resid = n - edf,
          change = change,
          R.sq = bR2,
          aic = aic
        )

        if (!is.null(keep)) {
          keep.list[[nm]] <- keep(fit, bR2)
        }
      },
      error = function(e) {
        if (trace > 0) {
          cat("Error updating model:", e$message, "\n")
        }
        break
      }
    )

    count.steps <- count.steps + 1
  }

  # Finalize the result
  if (!is.null(keep) && nm > 0) {
    fit$keep <- re.arrange(keep.list[seq(nm)])
  }

  fit <- make.step_gam(models = models[seq(nm)], fit, object)

  if (trace > 0) {
    print(fit$anova)
  }

  return(fit)
}
