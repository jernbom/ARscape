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
#'   Options: "genus", "species", "group".
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

  # Pre-allocate matrix for simulation results
  # Rows = representations, Cols = 1000 simulations
  sim_matrix <- matrix(NA, nrow = length(representations), ncol = 1000)
  rownames(sim_matrix) <- representations

  # 1. Generate Null Distributions
  for (i in seq_along(representations)) {
    n <- representations[i]

    # We replicate the mean calculation 1000 times
    sim_means <- replicate(1000, {
      mean(sample(pool_values, n, replace = FALSE))
    })

    sim_matrix[i, ] <- sim_means
  }

  # 2. Calculate Gaussian Parameters (Mean and SD)
  # For a normal distribution, the Maximum Likelihood Estimates (MLE) are exactly
  # the sample mean and standard deviation. No numerical fitting required!
  dist_info <- data.frame(
    total_peps = representations,
    mean = NA_real_,
    sd = NA_real_
  )
  fits <- list()

  for (i in seq_along(representations)) {
    valid_data <- sim_matrix[i, ]
    valid_data <- valid_data[!is.na(valid_data)]

    # Safety check for variance
    if (length(valid_data) > 1 && stats::var(valid_data) > 0) {
      dist_info$mean[i] <- mean(valid_data)
      dist_info$sd[i]   <- stats::sd(valid_data)
    } else {
      # Skip if there is no variance
      warning(paste("No variance in null data for n =", representations[i], "- using fallback values."))
      dist_info$mean[i] <- mean(valid_data)
      dist_info$sd[i]   <- 1e-6
    }

    # Store a pseudo-fit object to preserve compatibility with quantile extraction
    fits[[i]] <- list(
      data = valid_data,
      estimate = list(mean = dist_info$mean[i], sd = dist_info$sd[i])
    )
  }

  names(fits) <- representations
  dist_info$fits <- fits

  # 3. Spline Interpolation of Normal Parameters
  # Mean can be negative, so we interpolate on linear scale.
  # SD is strictly positive, so log10 remains safer to prevent predicting negative SDs.
  mean_spline <- stats::smooth.spline(x = log10(dist_info$total_peps), y = dist_info$mean, spar = 0.5)
  sd_spline   <- stats::smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$sd), spar = 0.5)

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
