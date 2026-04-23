#' Calculate Aggregate Reactivity Scores
#'
#' Generates aggregate reactivity scores by comparing the average fold change
#' of a group of peptides to null distributions derived from random sampling.
#' Assumes a Gaussian (Normal) null distribution.
#'
#' @param norm_log Data frame containing intermediate reactivity metrics.
#' @param all_peptide_fcs Data frame containing all peptide fold changes for an individual sample.
#' @param positives Data frame containing taxa determined to be significant (excluded from null).
#' @param exclusion_method Character string; method for excluding likely reactive peptides.
#'  Options: "genus", "species", "group".
#' @param seed Integer seed for reproducibility of random sampling.
#'
#' @return A data frame with added columns for mean, sd, ARscore, p_val, and nlog_p.
#'
#' @importFrom stats predict pnorm smooth.spline var sd dnorm
#' @importFrom dplyr filter mutate select
#' @export
calc_arscore <- function(norm_log,
                         all_peptide_fcs,
                         positives,
                         exclusion_method = "genus") {

  # Representations: number of peptides to simulate distributions for
  max_group_npep <- max(norm_log$total_peps)
  len_rep <- ceiling(log2(max_group_npep/30))
  representations <- 30*2^(0:len_rep)

  # Filter out known positives from the background pool
  if (exclusion_method %in% c("species", "group")) {
    background_pool <- all_peptide_fcs %>%
      dplyr::filter(!taxon_species %in% positives$taxon_species)
  } else {
    background_pool <- all_peptide_fcs %>%
      dplyr::filter(!taxon_genus %in% positives$taxon_genus)
  }

  pool_values <- background_pool$log2fc

  # Clean the source data
  if (any(!is.finite(pool_values))) {
    warning("NAs or Infinite values found in pool_values. Dropping prior to simulation.")
    pool_values <- pool_values[is.finite(pool_values)]
  }

  max_needed <- max(representations)
  if (length(pool_values) < max_needed) {
    warning(sprintf("Insufficient background peptides (Pool size: %d, Max needed: %d). Returning NAs.",
                    length(pool_values), max_needed))

    # Gracefully return NAs so the pipeline doesn't crash
    results <- norm_log %>%
      dplyr::mutate(
        null_mean = NA_real_,
        null_sd   = NA_real_,
        ARscore   = NA_real_,
        p_val     = NA_real_,
        nlog_p    = NA_real_
      )

    # Return empty dist_info to match expected list output
    return(list(results, data.frame()))
  }

  # Pre-allocate matrix for simulation results
  # Rows = representations, Cols = 1000 simulations
  sim_matrix <- matrix(NA, nrow = length(representations), ncol = 1000)
  rownames(sim_matrix) <- representations

  # 1. Generate Null Distributions
  for (i in seq_along(representations)) {
    n <- representations[i]

      # We replicate the mean calculation 1000 times
      sim_means <- replicate(1000, {
        # sample.int to prevent undesired behavior at the edge case length(pool_values) == 1
        mean(pool_values[base::sample.int(length(pool_values), n, replace = FALSE)])
      })

      sim_matrix[i, ] <- sim_means
  }

  # 2. Calculate Gaussian Parameters (Mean and SD)
  dist_info <- data.frame(
    total_peps = representations,
    mean = NA_real_,
    sd = NA_real_
  )
  fits <- list()

  for (i in seq_along(representations)) {
    valid_data <- sim_matrix[i, ]
    valid_data <- valid_data[!is.na(valid_data)]
    n_valid <- length(valid_data)

    # Logic tree:
    # 1. Compute mean and sd if possible (>1 valid points AND variance exists)
    if (n_valid > 1 && isTRUE(stats::sd(valid_data) > 1e-6) ) {
      dist_info$mean[i] <- mean(valid_data)
      dist_info$sd[i]   <- stats::sd(valid_data)

      # 2. Else mean only and supply small value for sd (exactly 1 point OR zero variance)
    } else if (n_valid >= 1) {
      warning(sprintf("No/undefined variance in null data for n = %d - using fallback SD.", representations[i]))
      dist_info$mean[i] <- mean(valid_data)
      dist_info$sd[i]   <- 1e-6

      # 3. Else skip the spline for that representation (0 valid points)
    } else {
      warning(sprintf("No valid simulation data for n = %d - skipping representation.", representations[i]))
      # Variables naturally remain NA_real_, which filters them out safely later
    }

    # Store a pseudo-fit object to preserve compatibility with quantile extraction
    fits[[i]] <- list(
      data = valid_data,
      estimate = list(mean = dist_info$mean[i], sd = dist_info$sd[i])
    )
  }

  names(fits) <- representations
  dist_info$fits <- fits

  # Sanitize dist_info before the spline
  dist_clean <- dist_info %>%
    dplyr::filter(is.finite(mean), is.finite(sd))

  # Enforce a minimum data point check
  # smooth.spline requires >= 4 unique x values
  if (nrow(dist_clean) < 4) {
    stop("Critical Error: Fewer than 4 valid data points remaining. Cannot fit smooth.spline(). Check background_pool for excessive NAs or zero variance.")
  }

  # 3. Spline Interpolation of Normal Parameters
  mean_spline <- stats::smooth.spline(x = log10(dist_clean$total_peps), y = dist_clean$mean, spar = 0.5)
  sd_spline   <- stats::smooth.spline(x = log10(dist_clean$total_peps), y = log10(dist_clean$sd), spar = 0.5)

  # 4. Apply to Data
  # Predict parameters for the actual 'total_peps' observed in the data
  results <- norm_log %>%
    dplyr::mutate(
      null_mean = stats::predict(mean_spline, log10(total_peps))$y,
      null_sd   = 10^(stats::predict(sd_spline, log10(total_peps))$y),

      # Z-score Calculation (Standard algebraic normal standardisation)
      ARscore = (score_norm - null_mean) / null_sd,

      # Calculate one-sided P-values
      p_val = stats::pnorm(ARscore, lower.tail = FALSE),
      nlog_p = -1 * stats::pnorm(ARscore, lower.tail = FALSE, log.p = TRUE) / log(10)
    )

  debug_results <- list(results, dist_info)

  return(debug_results)
}
