#' Calculate Z-scores from Generalized Gamma Distribution
#'
#' This function computes Z-scores (quantiles of the standard normal distribution)
#' corresponding to the cumulative probability of values in a Generalized Gamma distribution.
#' It uses the "split-tail" logic from `limma::zscoreGamma` to preserve numerical precision
#' for extreme values.
#'
#' @param q Numeric vector of quantiles (observed values).
#' @param mu Numeric vector of location parameters (from flexsurv::dgengamma).
#' @param sigma Numeric vector of scale parameters (must be positive).
#' @param Q Numeric vector of shape parameters.
#'
#' @return A numeric vector of Z-scores.
#' @keywords internal
#'
#' @importFrom stats qnorm
#' @importFrom flexsurv pgengamma qgengamma
zscoreGenGamma <- function(q, mu, sigma, Q) {
  # Initialize Z vector
  z <- q
  n <- length(q)

  # Recycle parameters to match length of q
  mu    <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  Q     <- rep_len(Q, n)

  # 1. Determine the "Pivot Point" (Median)
  # Splitting at the median ensures we always integrate the smaller tail,
  # preventing floating point underflow (1 - epsilon).
  median_val <- flexsurv::qgengamma(0.5, mu = mu, sigma = sigma, Q = Q)

  # 2. Identify Upper Tail values (Observed > Median)
  up <- (q > median_val)

  # 3. Handle Upper Tail (Right side)
  # Calculate P(X > q) using lower.tail = FALSE and log.p = TRUE
  if (any(up)) {
    log_p_upper <- flexsurv::pgengamma(q[up], mu = mu[up], sigma = sigma[up], Q = Q[up],
                                       lower.tail = FALSE, log.p = TRUE)

    # Convert log-probability directly to Z-score
    z[up] <- stats::qnorm(log_p_upper, lower.tail = FALSE, log.p = TRUE)
  }

  # 4. Handle Lower Tail (Left side)
  # Calculate P(X <= q) using lower.tail = TRUE and log.p = TRUE
  if (any(!up)) {
    log_p_lower <- flexsurv::pgengamma(q[!up], mu = mu[!up], sigma = sigma[!up], Q = Q[!up],
                                       lower.tail = TRUE, log.p = TRUE)

    z[!up] <- stats::qnorm(log_p_lower, lower.tail = TRUE, log.p = TRUE)
  }

  return(z)
}
