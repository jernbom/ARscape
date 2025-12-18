#' Calculate Augmented Reactivity Landscape Scores
#'
#' Generates aggregate reactivity scores by comparing the average fold change
#' of a group of peptides to null distributions derived from random sampling.
#'
#' @param norm_log Data frame containing intermediate reactivity metrics.
#' @param all_peptide_fcs Data frame containing all peptide fold changes for an individual sample.
#' @param positives Data frame containing taxa determined to be significant (excluded from null).
#' @param exclusion_method Character string; method for excluding likely reactive peptides.
#'   Options: "genus", "species", "group".
#' @param seed Integer seed for reproducibility of random sampling.
#'
#' @return A data frame with added columns for shape, rate, mean, variance, ARscore, and p_val.
#'
#' @importFrom stats lm predict pnorm smooth.spline var rnorm
#' @importFrom fitdistrplus fitdist
#' @importFrom dplyr filter mutate %>%
#' @importFrom limma zscoreGamma
#' @export
calculate_landscape <- function(norm_log,
                                all_peptide_fcs,
                                positives,
                                exclusion_method = "genus",
                                seed = 120) {

  # Set seed for reproducibility
  set.seed(seed)

  # Representations: number of peptides to simulate distributions for
  max_group_npep <- max(all_peptide_fcs$total_peps)
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

    # Handle absolute zeros to prevent log(0) issues in gamma fitting later
    sim_means[sim_means == 0] <- min(sim_means[sim_means != 0], na.rm = TRUE)

    sim_matrix[i, ] <- sim_means
  }

  # 2. Fit Gamma Distributions
  dist_info <- data.frame(
    total_peps = representations,
    shape = NA_real_,
    rate = NA_real_
  )

  for (i in seq_along(representations)) {
    valid_data <- sim_matrix[i, ]
    valid_data <- valid_data[valid_data > 0 & !is.na(valid_data)]

    # Safety check for variance
    if (length(valid_data) > 1 && stats::var(valid_data) > 0) {
      # Use tryCatch to prevent crashing on convergence failure
      fit <- tryCatch({
        fitdistrplus::fitdist(valid_data, "gamma", control = list(maxit = 1000))
      }, error = function(e) return(NULL))

      if (!is.null(fit)) {
        dist_info$shape[i] <- fit$estimate[["shape"]]
        dist_info$rate[i]  <- fit$estimate[["rate"]]
      } else {
        warning(paste("Gamma fit failed for n =", representations[i]))
        dist_info$shape[i] <- 1
        dist_info$rate[i] <- 1
      }
    } else {
      dist_info$shape[i] <- 1
      dist_info$rate[i] <- 1
    }
  }

  # 3. Spline Interpolation (Smoothing parameters across N)
  # We work in log10 space for linearity
  log_peps <- log10(dist_info$total_peps)

  # Fit splines
  shape_spline <- stats::smooth.spline(x = log_peps, y = log10(dist_info$shape), spar = 0.5)
  rate_spline  <- stats::smooth.spline(x = log_peps, y = log10(dist_info$rate), spar = 0.5)

  # 4. Apply to Data
  # Predict parameters for the actual 'total_peps' observed in the data
  results <- norm_log %>%
    dplyr::mutate(
      pred_log_peps = log10(total_peps),
      shape_pred = 10^(stats::predict(shape_spline, pred_log_peps)$y),
      rate_pred  = 10^(stats::predict(rate_spline, pred_log_peps)$y),

      # Metrics
      mean_expected = shape_pred / rate_pred,

      # Z-score Calculation
      ARscore = limma::zscoreGamma(score_norm, shape = shape_pred, rate = rate_pred),
      # Calculate P-values
      p_val = stats::pnorm(ARscore, lower.tail = FALSE),
      nlog_p = -log10(p_val)
    )

  return(results)
}

