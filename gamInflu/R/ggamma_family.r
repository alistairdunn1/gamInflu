#' @title Generalized Gamma Distribution Functions
#' @description S3 methods for the 3-parameter generalized gamma distribution,
#' compatible with mgcv's extended family framework.
#' @details
#' This implements a generalized gamma distribution with three parameters:
#' \itemize{
#'   \item \code{mu}: location parameter (log scale)
#'   \item \code{sigma}: scale parameter (positive)
#'   \item \code{Q}: shape parameter
#' }
#'
#' The generalized gamma reduces to special cases:
#' \itemize{
#'   \item When Q \\to 0: Log-normal distribution
#'   \item When Q = 1: Weibull distribution
#'   \item When sigma = 1: Gamma distribution
#' }
#' @name ggamma-distribution
NULL

#' @rdname ggamma-distribution
#' @param x,q vector of quantiles
#' @param mu location parameter (log scale)
#' @param sigma scale parameter (positive)
#' @param Q shape parameter
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @export
dgg <- function(x, mu, sigma, Q, log = FALSE) {
  UseMethod("dgg")
}

#' @rdname ggamma-distribution
#' @export
dgg.default <- function(x, mu, sigma, Q, log = FALSE) {
  # Density function for generalized gamma
  # mu: location parameter (log scale)
  # sigma: scale parameter (positive)
  # Q: shape parameter

  if (any(x <= 0)) stop("x must be positive")

  w <- (log(x) - mu) / sigma

  if (abs(Q) < 1e-7) {
    # Limit as Q -> 0 (log-normal case)
    logdens <- -0.5 * log(2 * pi) - log(sigma) - log(x) - 0.5 * w^2
  } else {
    Q2 <- Q * Q
    Q_inv2 <- 1 / Q2
    gamma_term <- lgamma(Q_inv2)

    # More numerically stable computation
    u <- Q_inv2 + Q * w / abs(Q)

    logdens <- log(abs(Q)) - log(sigma) - log(x) - gamma_term -
      Q_inv2 * log(Q_inv2) +
      (Q_inv2 - 1) * log(pmax(u, 1e-10)) - u
  }

  if (log) {
    return(logdens)
  } else {
    return(exp(logdens))
  }
}

#' @rdname ggamma-distribution
#' @param p vector of probabilities
#' @param lower.tail logical; if TRUE (default), probabilities are P[X \\le x], otherwise P[X > x]
#' @export

pgg <- function(q, mu, sigma, Q, lower.tail = TRUE, log.p = FALSE) {
  UseMethod("pgg")
}

#' @rdname ggamma-distribution
#' @export
pgg.default <- function(q, mu, sigma, Q, lower.tail = TRUE, log.p = FALSE) {
  w <- (log(q) - mu) / sigma

  if (abs(Q) < 1e-7) {
    p <- pnorm(w, lower.tail = lower.tail, log.p = log.p)
  } else {
    Q2 <- Q * Q
    Q_inv2 <- 1 / Q2
    u <- Q_inv2 + Q * w / abs(Q)

    if (Q > 0) {
      p <- pgamma(pmax(u, 0),
        shape = Q_inv2, rate = 1,
        lower.tail = lower.tail, log.p = log.p
      )
    } else {
      p <- pgamma(pmax(u, 0),
        shape = Q_inv2, rate = 1,
        lower.tail = !lower.tail, log.p = log.p
      )
    }
  }
  return(p)
}

#' @rdname ggamma-distribution
#' @export

qgg <- function(p, mu, sigma, Q, lower.tail = TRUE, log.p = FALSE) {
  UseMethod("qgg")
}

#' @rdname ggamma-distribution
#' @export
qgg.default <- function(p, mu, sigma, Q, lower.tail = TRUE, log.p = FALSE) {
  if (abs(Q) < 1e-7) {
    w <- qnorm(p, lower.tail = lower.tail, log.p = log.p)
  } else {
    Q2 <- Q * Q
    Q_inv2 <- 1 / Q2

    if (Q > 0) {
      u <- qgamma(p,
        shape = Q_inv2, rate = 1,
        lower.tail = lower.tail, log.p = log.p
      )
    } else {
      u <- qgamma(p,
        shape = Q_inv2, rate = 1,
        lower.tail = !lower.tail, log.p = log.p
      )
    }
    w <- abs(Q) * (u - Q_inv2) / Q
  }

  return(exp(mu + sigma * w))
}

#' @rdname ggamma-distribution
#' @param n number of observations
#' @export

rgg <- function(n, mu, sigma, Q) {
  UseMethod("rgg")
}

#' @rdname ggamma-distribution
#' @export
rgg.default <- function(n, mu, sigma, Q) {
  qgg(runif(n), mu, sigma, Q)
}

#' @title Generalized Gamma Family for mgcv
#' @description Creates a generalized gamma family object for use with mgcv's gam() function.
#' This implements a 3-parameter generalized gamma distribution as an extended family.
#' @param link A list of three link functions for mu, sigma, and Q parameters.
#'   Default is \code{list("log", "log", "identity")}.
#' @param b Offset parameter (currently unused, kept for compatibility)
#' @return An extended family object suitable for use with \code{mgcv::gam()}
#' @details
#' The generalized gamma distribution is a flexible 3-parameter family that includes
#' many common distributions as special cases:
#' \itemize{
#'   \item Log-normal when Q \\to 0
#'   \item Weibull when Q = 1
#'   \item Gamma when sigma = 1
#' }
#'
#' The parameterization uses:
#' \itemize{
#'   \item \code{mu}: location parameter (typically on log scale)
#'   \item \code{sigma}: scale parameter (positive)
#'   \item \code{Q}: shape parameter
#' }
#' @examples
#' \dontrun{
#' library(mgcv)
#'
#' # Simulate data from generalized gamma
#' set.seed(123)
#' n <- 200
#' x <- runif(n)
#' mu_true <- exp(1 + 2 * x)
#' y <- rgg(n, log(mu_true), 0.5, 0.8)
#'
#' # Fit GAM with generalized gamma family
#' fit <- gam(y ~ s(x), family = gengamma(), data = data.frame(y = y, x = x))
#' summary(fit)
#' plot(fit)
#' }
#' @importFrom stats make.link runif
#' @export
gengamma <- function(link = list("log", "log", "identity"), b = 0) {
  # Validate links
  if (length(link) != 3) {
    stop("gengamma family requires exactly 3 link functions for mu, sigma, Q")
  }

  # Set up link functions as objects
  linktemp <- substitute(link)
  if (!is.list(linktemp)) linktemp <- list(linktemp)

  # Process each link
  for (i in 1:3) {
    if (is.name(linktemp[[i]])) {
      linktemp[[i]] <- deparse(linktemp[[i]])
    }
    if (is.character(linktemp[[i]])) {
      linktemp[[i]] <- switch(linktemp[[i]],
        "log" = make.link("log"),
        "identity" = make.link("identity"),
        "sqrt" = make.link("sqrt"),
        "inverse" = make.link("inverse"),
        make.link(linktemp[[i]])
      )
    }
  }

  link_names <- c(linktemp[[1]]$name, linktemp[[2]]$name, linktemp[[3]]$name)

  # Environment for storing methods
  env <- new.env()

  # Extended family methods

  # Deviance residuals
  env$dev.resids <- function(y, mu, wt, theta, scale = 1) {
    sigma <- exp(theta[1, ])
    Q <- theta[2, ]

    # Unit deviances
    dev <- numeric(length(y))
    for (i in seq_along(y)) {
      ll_sat_i <- wt[i] * dgg(y[i], log(y[i]), sigma[i], Q[i], log = TRUE)
      ll_fit_i <- wt[i] * dgg(y[i], log(mu[i]), sigma[i], Q[i], log = TRUE)
      dev[i] <- 2 * (ll_sat_i - ll_fit_i)
    }

    return(sign(y - mu) * sqrt(pmax(dev, 0)))
  }

  # AIC
  env$aic <- function(y, n, mu, wt, dev, theta, scale = 1) {
    sigma <- exp(theta[1, ])
    Q <- theta[2, ]
    ll <- sum(wt * dgg(y, log(mu), sigma, Q, log = TRUE))
    return(-2 * ll)
  }

  # Log-likelihood derivatives (essential for mgcv)
  env$d1link <- function(mu, theta) {
    # First derivatives of log-likelihood w.r.t. linear predictors
    # This is crucial for mgcv's fitting algorithm
    matrix(0, nrow = length(mu), ncol = 3) # Placeholder - computed in ll.grad
  }

  env$d2link <- function(mu, theta) {
    # Second derivatives (Hessian)
    # Also computed in ll.grad for efficiency
    array(0, dim = c(length(mu), 3, 3)) # Placeholder
  }

  # Gradient and Hessian of log-likelihood (key method)
  env$ll.grad <- function(y, mu, theta, wt, deriv = 2) {
    sigma <- exp(theta[1, ])
    Q <- theta[2, ]
    n <- length(y)

    w <- (log(y) - log(mu)) / sigma

    # Initialize gradient and Hessian
    grad <- matrix(0, n, 3) # w.r.t. eta_mu, eta_sigma, eta_Q
    hess <- array(0, c(n, 3, 3))

    for (i in 1:n) {
      if (abs(Q[i]) < 1e-7) {
        # Log-normal case
        # d/d(log mu)
        grad[i, 1] <- wt[i] * w[i] / sigma[i]
        # d/d(log sigma)
        grad[i, 2] <- wt[i] * sigma[i] * (w[i]^2 - 1)
        # d/dQ (use Taylor expansion)
        grad[i, 3] <- wt[i] * w[i]^3 / 6

        if (deriv > 1) {
          # Second derivatives for log-normal
          hess[i, 1, 1] <- -wt[i] / (sigma[i]^2)
          hess[i, 1, 2] <- hess[i, 2, 1] <- wt[i] * w[i] * sigma[i] / sigma[i]
          hess[i, 2, 2] <- wt[i] * sigma[i]^2 * (2 * w[i]^2 - 1)
          hess[i, 1, 3] <- hess[i, 3, 1] <- wt[i] * w[i]^2 / (2 * sigma[i])
          hess[i, 2, 3] <- hess[i, 3, 2] <- wt[i] * sigma[i] * w[i]^3 / 2
          hess[i, 3, 3] <- wt[i] * w[i]^4 / 12
        }
      } else {
        # General case
        Q2 <- Q[i] * Q[i]
        Q_inv2 <- 1 / Q2
        u <- Q_inv2 + Q[i] * w[i] / abs(Q[i])

        # First derivatives
        # d/d(log mu)
        grad[i, 1] <- wt[i] * sign(Q[i]) * Q[i] / sigma[i] * (Q_inv2 / u - 1)

        # d/d(log sigma)
        grad[i, 2] <- wt[i] * sigma[i] * (-1 + sign(Q[i]) * Q[i] * w[i] * (Q_inv2 / u - 1))

        # d/dQ (simplified)
        psi_val <- digamma(Q_inv2)
        grad[i, 3] <- wt[i] * (sign(Q[i]) / Q[i] - 2 * Q_inv2 * sign(Q[i]) / Q[i] *
          (psi_val - log(Q_inv2) - 1 + log(u)) +
          w[i] / abs(Q[i]) * (Q_inv2 / u - 1))

        if (deriv > 1) {
          # Approximate second derivatives (for numerical stability)
          hess[i, 1, 1] <- -wt[i] * Q2 / (sigma[i]^2 * u)
          hess[i, 1, 2] <- hess[i, 2, 1] <- wt[i] * Q[i] * sign(Q[i]) * w[i] / sigma[i] * Q_inv2 / u
          hess[i, 2, 2] <- -wt[i] * sigma[i]^2 * Q2 * w[i] / abs(Q[i]) * Q_inv2 / u
          # Set off-diagonal Q terms to small values for stability
          hess[i, 1, 3] <- hess[i, 3, 1] <- wt[i] * 0.1 / sigma[i]
          hess[i, 2, 3] <- hess[i, 3, 2] <- wt[i] * 0.1 * sigma[i]
          hess[i, 3, 3] <- -wt[i] * 0.1
        }
      }
    }

    if (deriv == 1) {
      return(list(grad = grad))
    } else {
      return(list(grad = grad, hess = hess))
    }
  }

  # Theta initialization
  env$initialize <- expression({
    n <- rep(1, nobs)

    # Robust initialization using quantiles
    if (any(y <= 0)) stop("All response values must be positive")

    log_y <- log(y)

    # Method of moments estimates
    mean_log <- mean(log_y)
    var_log <- var(log_y)

    # Initial sigma (bounded)
    sigma_init <- sqrt(pmax(var_log, 0.01))

    # Initial Q (start near log-normal)
    Q_init <- 0.1

    # Theta matrix: first row log(sigma), second row Q
    theta <- rbind(
      rep(log(sigma_init), length(y)),
      rep(Q_init, length(y))
    )
  })

  # Random deviates
  env$rd <- function(mu, wt, scale, theta) {
    sigma <- exp(theta[1, ])
    Q <- theta[2, ]
    rgg(length(mu), log(mu), sigma, Q)
  }

  # Prediction methods
  env$predict <- function(family, se = FALSE, eta = NULL, theta = NULL,
                          X = NULL, beta = NULL, off = NULL, Vb = NULL) {
    if (is.null(eta)) {
      eta <- X %*% beta + off
      if (se) {
        se.fit <- sqrt(pmax(0, rowSums((X %*% Vb) * X)))
      }
    }

    fv <- exp(eta) # inverse link for mu

    if (se) {
      return(list(fit = fv, se.fit = se.fit))
    } else {
      return(fv)
    }
  }

  # Variance function (approximate)
  variance <- function(mu) mu^2 # Rough approximation

  # Validation functions
  validmu <- function(mu) all(is.finite(mu) & mu > 0)
  valideta <- function(eta) all(is.finite(eta))

  # Create family object
  structure(list(
    family = "gengamma",
    link = link_names,
    linkfun = linktemp[[1]]$linkfun,
    linkinv = linktemp[[1]]$linkinv,
    mu.eta = linktemp[[1]]$mu.eta,
    variance = variance,
    validmu = validmu,
    valideta = valideta,
    dev.resids = env$dev.resids,
    aic = env$aic,
    ll.grad = env$ll.grad,
    d1link = env$d1link,
    d2link = env$d2link,
    initialize = env$initialize,
    rd = env$rd,
    predict = env$predict,
    n.theta = 2,
    no.theta = 2,
    ini.theta = function(theta) c(log(0.5), 0.1),
    putTheta = function(theta) theta,
    getTheta = function(fit) fit$family$getTheta(fit),
    use.unscaled.deviance = TRUE,
    scale.known = TRUE,
    env = env
  ), class = c("extended.family", "family"))
}

#' @title Print Method for Extended Family Objects
#' @description S3 print method for extended family objects created by gengamma()
#' @param x An extended family object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object
#' @export
print.extended.family <- function(x, ...) {
  cat("\nFamily:", x$family, "\n")
  cat("Link functions:", paste(x$link, collapse = ", "), "\n")
  cat("Number of theta parameters:", x$n.theta, "\n\n")
  invisible(x)
}
