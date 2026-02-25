#' Calculate Aggregate Reactivity Scores
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

    # Handle absolute zeros to prevent log(0) issues in gamma fitting later
    sim_means[sim_means == 0] <- min(sim_means[sim_means != 0], na.rm = TRUE)/2

    sim_matrix[i, ] <- sim_means
  }

  # 2. Fit Gamma Distributions
  dist_info <- data.frame(
    total_peps = representations,
    shape = NA_real_,
    rate = NA_real_
  )
  gammafits <- list()

  for (i in seq_along(representations)) {
    valid_data <- sim_matrix[i, ]
    valid_data <- valid_data[valid_data > 0 & !is.na(valid_data)]

    # Safety check for variance
    if (length(valid_data) > 1 && stats::var(valid_data) > 0) {

      # Use tryCatch for errors and withCallingHandlers for targeted warning suppression
      fit <- tryCatch({
        withCallingHandlers(
          expr = fitdistrplus::fitdist(valid_data, "gamma", control = list(maxit = 1000)),
          warning = function(w) {
            # Silently suppress the NaN warning and the cov2cor warning
            if (grepl("NaNs produced|diag\\(V\\) had non-positive", conditionMessage(w))) {
              invokeRestart("muffleWarning")
            }
            # All other warnings bypass this and print natively
          }
        )
      }, error = function(e) {
        # Errors print to console and break execution
        message(sprintf("Error during Gamma fitting for n = %d: %s", representations[i], conditionMessage(e)))
        stop(e)
        # NOTE: If you prefer to preserve Rule 5 (the fallback to shape/rate = 1) instead
        # of crashing the pipeline, comment out `stop(e)` above and uncomment `return(NULL)` below:
        # return(NULL)
      })

      if (!is.null(fit)) {
        dist_info$shape[i] <- fit$estimate[["shape"]]
        dist_info$rate[i]  <- fit$estimate[["rate"]]
        gammafits[[i]] <- fit
      } else {
        # 5. Warning issued if shape and rate are Null (Fallback logic, currently inactive as tryCatch currently stops on error, see above.)
        warning(paste("Gamma fit failed for n =", representations[i]))
        dist_info$shape[i] <- 1
        dist_info$rate[i] <- 1
        gammafits[[i]] <- list()
      }
    } else {
      # 6. Skip fitting if there is no variance
      warning(paste("Skipping Gamma fitting for n =", representations[i], "- No variance in data."))
      dist_info$shape[i] <- 1
      dist_info$rate[i] <- 1
      gammafits[[i]] <- list()
    }
  }

  names(gammafits) <- representations

  dist_info <- dist_info %>% mutate(mean = shape / rate) %>%
    mutate(variance = shape / rate^2)

  dist_info[["gammafits"]] <- gammafits

  # 3. Spline Interpolation of Gamma Parameters
  shape_spline <- stats::smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$shape), spar = 0.5)
  rate_spline  <- stats::smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$rate), spar = 0.5)

  # 4. Apply to Data
  # Predict parameters for the actual 'total_peps' observed in the data
  results <- norm_log %>%
    dplyr::mutate(
      shape = 10^(stats::predict(shape_spline, log10(total_peps))$y),
      rate  = 10^(stats::predict(rate_spline, log10(total_peps))$y),

      # Z-score Calculation
      ARscore = limma::zscoreGamma(score_norm, shape = shape, rate = rate),
      # Calculate one-sided P-values
      p_val = stats::pnorm(ARscore, lower.tail = FALSE),
      nlog_p = -1 * stats::pnorm(ARscore, lower.tail = FALSE, log.p = TRUE) / log(10)
    )

  debug_results <- list(results, dist_info)

  return(debug_results)
}
