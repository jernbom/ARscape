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
#' @return A data frame with added columns for mu, sigma, Q, ARscore, and p_val.
#'
#' @importFrom stats lm predict pnorm smooth.spline var sd quantile
#' @importFrom fitdistrplus fitdist
#' @importFrom flexsurv pgengamma dgengamma
#' @importFrom dplyr filter mutate %>% select
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
  sim_matrix <- matrix(NA, nrow = length(representations), ncol = 1000)
  rownames(sim_matrix) <- representations

  # 1. Generate Null Distributions
  for (i in seq_along(representations)) {
    n <- representations[i]
    sim_means <- replicate(1000, {
      mean(sample(pool_values, n, replace = FALSE))
    })

    # Handle absolute zeros/negatives if necessary (GenGamma requires x > 0)
    # Assuming positive data or applying a shift/filter as per original code
    sim_means[sim_means <= 0] <- min(sim_means[sim_means > 0], na.rm = TRUE)/2
    sim_matrix[i, ] <- sim_means
  }

  # 2. Fit Generalized Gamma Distributions
  dist_info <- data.frame(
    total_peps = representations,
    mu = NA_real_,
    sigma = NA_real_,
    Q = NA_real_
  )
  gammafits <- list()

  for (i in seq_along(representations)) {
    # 1. Prepare Data
    valid_data <- sim_matrix[i, ]
    valid_data <- valid_data[valid_data > 0 & !is.na(valid_data)]

    # Calculate Diagnostics immediately for reporting
    n_pts   <- length(valid_data)
    var_val <- if(n_pts > 1) stats::var(valid_data) else NA
    mean_val <- if(n_pts > 0) mean(valid_data) else NA
    cv_val  <- if(!is.na(var_val) && mean_val != 0) sqrt(var_val)/mean_val else NA

    fit <- NULL
    fail_reason <- NULL

    # 2. Check Prerequisites (Basic Length Check only)
    if (length(valid_data) <= 5) {
      fail_reason <- "Insufficient data length (<= 5 points)."
    } else {

      # 3. Attempt Fit with Warning/Error Capture
      fit <- tryCatch({
        withCallingHandlers({

          # Initialize with Log-Normal moments (Q=0 approximation)
          # GenGamma is sensitive; good start values are critical.
          start_vals <- list(mu = mean(log(valid_data)),
                             sigma = sd(log(valid_data)),
                             Q = 0)

          # Fit using L-BFGS-B to constrain Sigma > 0
          fitdistrplus::fitdist(valid_data, "gengamma",
                                start = start_vals,
                                optim.method = "L-BFGS-B",
                                lower = c(-Inf, 1e-5, -Inf), # mu, sigma, Q
                                control = list(maxit = 2000))

        }, warning = function(w) {
          # --- WARNING FILTER ---
          call_context <- if (!is.null(w$call)) paste(deparse(w$call), collapse = " ") else "Unknown"

          # We log helpful stats even in warnings now
          warn_msg <- paste0(
            "[WARNING] Sample: ", sample_id, " | N-peps: ", representations[i],
            "\n   Msg: ", w$message,
            "\n   Stats -> Var: ", format(var_val, digits=4),
            " | CV: ", format(cv_val, digits=4),
            " | Mean: ", format(mean_val, digits=4)
          )

          # Muffle standard numerical warnings if they aren't fatal
          if (grepl("NaNs produced", w$message) || grepl("Hessian", w$message)) {
            invokeRestart("muffleWarning")
          } else {
            message(warn_msg)
            invokeRestart("muffleWarning")
          }
        })

      }, error = function(e) {
        return(e$message) # Return error string to be handled below
      })

      if (is.character(fit)) {
        fail_reason <- fit
        fit <- NULL
      }
    }

    # 4. Process Outcomes
    if (!is.null(fit)) {
      # --- SUCCESS ---
      dist_info$mu[i]    <- fit$estimate[["mu"]]
      dist_info$sigma[i] <- fit$estimate[["sigma"]]
      dist_info$Q[i]     <- fit$estimate[["Q"]]
      gammafits[[i]] <- fit

    } else {
      # --- FAILURE ---
      # Generate detailed debug report
      debug_text <- utils::capture.output({
        cat("\n--- GENGAMMA FIT FAILURE ---\n")
        cat("Reason:", fail_reason, "\n")
        cat("Sample ID:", sample_id, "| N-peptides:", representations[i], "\n")
        cat("Stats -> Mean:", mean_val, "| Var:", var_val, "| CV:", cv_val, "\n")
        cat("Summary of valid_data:\n")
        print(summary(valid_data))
        cat("----------------------------\n")
      })

      message(paste(debug_text, collapse = "\n"))

      # Leave parameters as NA so Spline skips them
      gammafits[[i]] <- list()
    }
  }

  names(gammafits) <- representations
  dist_info[["gammafits"]] <- gammafits

  # 3. Spline Interpolation
  # Filter to only successful fits for spline creation
  valid_fits <- dist_info %>% dplyr::filter(!is.na(mu))

  if (nrow(valid_fits) < 2) {
    stop("Not enough successful fits to generate splines. Check data variance.")
  }

  # Mu: Linear Spline
  mu_spline <- stats::smooth.spline(x = log10(valid_fits$total_peps),
                                    y = valid_fits$mu, spar = 0.5)

  # Sigma: Log Spline (Enforce positivity)
  sigma_spline <- stats::smooth.spline(x = log10(valid_fits$total_peps),
                                       y = log10(valid_fits$sigma), spar = 0.5)

  # Q: Linear Spline (Allow negative values)
  Q_spline <- stats::smooth.spline(x = log10(valid_fits$total_peps),
                                   y = valid_fits$Q, spar = 0.5)

  # 4. Apply to Data
  results <- norm_log %>%
    dplyr::mutate(
      # Predict parameters
      mu_pred    = stats::predict(mu_spline, log10(total_peps))$y,
      sigma_pred = 10^(stats::predict(sigma_spline, log10(total_peps))$y), # Undo log10
      Q_pred     = stats::predict(Q_spline, log10(total_peps))$y,

      # Z-score Calculation (Using your custom function)
      ARscore = zscoreGenGamma(score_norm, mu = mu_pred, sigma = sigma_pred, Q = Q_pred),

      # Calculate P-values (External calculation via pnorm is safer/cleaner)
      p_val = stats::pnorm(ARscore, lower.tail = FALSE),

      # Robust Negative Log P-value
      nlog_p = -1 * stats::pnorm(ARscore, lower.tail = FALSE, log.p = TRUE) / log(10)
    )

  debug_results <- list(results, dist_info, sim_matrix)

  return(debug_results)
}
