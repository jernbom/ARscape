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
#' @param tau numeric vector of hidden mean parameters (default 0 for standard SN)
#' @return numeric vector of standard normal Z-scores
zscoreSN <- function(q, xi = 0, omega = 1, alpha = 0, tau = 0) {

  # Initialize Z-score vector
  z <- q
  n <- length(q)

  # Recycle parameters to match the length of the input quantiles
  xi    <- rep_len(xi, length.out = n)
  omega <- rep_len(omega, length.out = n)
  alpha <- rep_len(alpha, length.out = n)
  tau   <- rep_len(tau, length.out = n)

  # Calculate the theoretical mean of the extended skew-normal distribution
  # Formula: mu = xi + omega * delta * (phi(tau_0) / Phi(tau_0))
  # where delta = alpha / sqrt(1 + alpha^2) and tau_0 = tau / sqrt(1 + alpha^2)
  # Note: When tau = 0, phi(0)/Phi(0) simplifies exactly to sqrt(2/pi)
  delta <- alpha / sqrt(1 + alpha^2)
  tau_0 <- tau / sqrt(1 + alpha^2)
  mean_sn <- xi + omega * delta * (dnorm(tau_0) / pnorm(tau_0))

  # Split the data into upper and lower halves based on the mean
  up <- (q > mean_sn)

  # --- LOWER TAIL ---
  # For values below the mean, standard CDF works perfectly.
  if (any(!up)) {
    p_lower <- sn::psn(q[!up], xi = xi[!up], omega = omega[!up], alpha = alpha[!up], tau = tau[!up])

    # Convert probability to Z-score
    z[!up] <- qnorm(p_lower, lower.tail = TRUE)
  }

  # --- UPPER TAIL ---
  # For values above the mean, we avoid `1 - psn(...)` to prevent floating-point
  # cancellation (which would artificially cap max Z-scores to ~8).
  # Instead, we use the symmetry property: P(X > x) == P(-X < -x)
  # where -X follows a Skew-Normal with negated location and slant (tau remains unchanged).
  if (any(up)) {
    p_upper <- sn::psn(-q[up], xi = -xi[up], omega = omega[up], alpha = -alpha[up], tau = tau[up])

    # Convert the upper-tail probability directly to a positive Z-score
    z[up] <- qnorm(p_upper, lower.tail = FALSE)
  }

  return(z)
}
