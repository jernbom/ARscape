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
                         sample_id = "Unknown",
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
    # 1. Prepare Data
    valid_data <- sim_matrix[i, ]
    valid_data <- valid_data[valid_data > 0 & !is.na(valid_data)]

    # Placeholders for logic flow
    fit <- NULL
    fail_reason <- NULL

    # 2. Check Prerequisites
    if (length(valid_data) <= 1 || stats::var(valid_data) <= 0) {
      fail_reason <- "Skipping Gamma fitting - No variance in data or insufficient length."
    } else {

      # 3. Attempt Fit with Warning Capture
      fit <- tryCatch({
        withCallingHandlers({
          # Fit the model
          fitdistrplus::fitdist(valid_data, "gamma",
                                start = list(shape = 1, rate = 1),
                                control = list(maxit = 1000))
        }, warning = function(w) {
          # --- WARNING FILTER ---
          call_context <- if (!is.null(w$call)) paste(deparse(w$call), collapse = " ") else "Unknown"
          is_nan_msg   <- w$message == "NaNs produced"
          is_dist_func <- grepl("dgamma", call_context) || grepl("pgamma", call_context)

          if (is_nan_msg && is_dist_func) {
            # SILENT PATH: Muffle the noise and continue
            invokeRestart("muffleWarning")
          } else {
            # LOUD PATH: This is an unexpected warning, so we log it!
            warn_msg <- paste0(
              "[WARNING] Sample: ", sample_id,
              " | N-peps: ", representations[i],
              "\n   Context: In ", call_context,
              "\n   Msg: ", w$message,
            )
            message(warn_msg)
            invokeRestart("muffleWarning")
          }
        })

      }, error = function(e) {
        # If it crashes completely, return the message string
        return(e$message)
      })

      # If 'fit' is a character string, it means tryCatch caught an error
      if (is.character(fit)) {
        fail_reason <- fit
        fit <- NULL
      }
    }

    # 4. Process Outcomes
    if (!is.null(fit)) {
      # --- SUCCESS ---
      dist_info$shape[i] <- fit$estimate[["shape"]]
      dist_info$rate[i]  <- fit$estimate[["rate"]]
      gammafits[[i]] <- fit

    } else {
      # --- FAILURE (From either Variance check OR fitdist crash) ---

      # A. Capture the Debug Info
      debug_text <- utils::capture.output({
        cat("\n--- GAMMA FIT FAILURE ---\n")
        cat("Reason:", fail_reason, "\n")
        cat("Sample ID:", sample_id, "| N-peptides:", representations[i], "\n")

        cat("Length of valid data:", length(valid_data), "\n")
        cat("Variance of valid data:", if(length(valid_data) > 1) stats::var(valid_data) else "NA", "\n")

        cat("Head of filtered data:\n")
        print(head(valid_data))

        cat("Summary of full simulation row (sim_matrix):\n")
        # coerce to numeric to ensure summary works nicely
        print(summary(as.numeric(sim_matrix[i, ])))
        cat("---------------------------\n")
      })

      # B. Send to Main Log (using message to be parallel-safe)
      message(paste(debug_text, collapse = "\n"))

      # C. Set Fallbacks
      dist_info$shape[i] <- 1
      dist_info$rate[i]  <- 1
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
      nlog_p = -log10(p_val)
    )

  debug_results <- list(results, dist_info, sim_matrix)

  return(debug_results)
}

