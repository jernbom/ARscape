#' Compute Z-scores based on Skew-Normal distributions
#'
#' Analogous to limma::zscoreGamma, this function converts quantiles from
#' a skew-normal distribution into standard normal Z-scores. It splits the
#' calculations at the mean to ensure precision at extreme values in both tails.
#'
#' @param q numeric vector of quantiles (e.g., your observed values)
#' @param xi numeric vector of location parameters
#' @param omega numeric vector of scale parameters
#' @param alpha numeric vector of slant (skew) parameters
#' @return numeric vector of standard normal Z-scores
zscoreSN <- function(q, xi = 0, omega = 1, alpha = 0) {

  # Initialize Z-score vector
  z <- q
  n <- length(q)

  # Recycle parameters to match the length of the input quantiles
  xi    <- rep_len(xi, length.out = n)
  omega <- rep_len(omega, length.out = n)
  alpha <- rep_len(alpha, length.out = n)

  # Calculate the theoretical mean of the skew-normal distribution for each observation
  # Formula: mu = xi + omega * delta * sqrt(2/pi)
  # where delta = alpha / sqrt(1 + alpha^2)
  delta <- alpha / sqrt(1 + alpha^2)
  mean_sn <- xi + omega * delta * sqrt(2 / pi)

  # Split the data into upper and lower halves based on the mean
  up <- (q > mean_sn)

  # --- LOWER TAIL ---
  # For values below the mean, standard CDF works perfectly.
  if (any(!up)) {
    p_lower <- sn::psn(q[!up], xi = xi[!up], omega = omega[!up], alpha = alpha[!up])

    # Convert probability to Z-score
    z[!up] <- qnorm(p_lower, lower.tail = TRUE)
  }

  # --- UPPER TAIL ---
  # For values above the mean, we avoid `1 - psn(...)` to prevent floating-point
  # cancellation (which would artificially cap max Z-scores to ~8).
  # Instead, we use the symmetry property: P(X > x) == P(-X < -x)
  # where -X follows a Skew-Normal with negated location and slant.
  if (any(up)) {
    p_upper <- sn::psn(-q[up], xi = -xi[up], omega = omega[up], alpha = -alpha[up])

    # Convert the upper-tail probability directly to a positive Z-score
    z[up] <- qnorm(p_upper, lower.tail = FALSE)
  }

  return(z)
}
