#' Run ARscape Analysis Pipeline
#'
#' The main wrapper function to process PhIP-Seq fold-change data and generate
#' Aggregate Reactivity Scores (ARscores).
#'
#' @section Parallelization & Progress:
#' This function leverages the `furrr` package for parallel processing and `progressr`
#' for progress reporting.
#'
#' To enable parallel execution:
#' \preformatted{
#' future::plan(future::multisession, workers = parallelly::availableCores())
#' }
#'
#' To enable progress reporting:
#' \preformatted{
#' progressr::with_progress({
#'   run_arscape(...)
#' })
#' }
#'
#' @param fold_change Data frame containing peptide-level fold changes.
#'   Must contain annotation columns (specified in `annotation_cols`) and sample columns.
#' @param beads_controls Optional data frame of beads-only controls. If provided,
#'   peptides with fold change > 1 in these controls are excluded.
#'   Alternatively, provide a character vector of `bad_peptides` to skip this calculation.
#' @param bad_peptides Optional character vector of peptide sequences/IDs to exclude manually.
#' @param annotation_cols Character vector of column names that are NOT samples.
#'   Defaults to standard PhIP-seq annotations: `c("pep_aa", "taxon_species", "taxon_genus")`.
#' @param max_iterations Integer; maximum number of iterations for the background definition algorithm.
#' @param p_cutoff Numeric; -log10 p-value cutoff for defining positive hits. Default is 4 (-log10(0.0001)).
#' @param score_cutoff Numeric; ARscore cutoff for defining positive hits. Default is 2.
#' @param min_peptides Integer; minimum number of peptides required per group to calculate a score. Default is 50.
#' @param exclusion_method Character; method for excluding reactive peptides during iteration.
#'   Options: "genus", "species", "group".
#'
#' @return A data frame containing the final ARscores and p-values for all samples.
#'
#' @importFrom dplyr select filter group_by summarise mutate left_join anti_join pick everything all_of any_of n distinct bind_rows pull
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map_dfr
#' @importFrom furrr future_map_dfr furrr_options
#' @importFrom progressr progressor
#' @export
run_arscape <- function(fold_change,
                        beads_controls = NULL,
                        bad_peptides = NULL,
                        annotation_cols = c("pep_aa", "taxon_species", "taxon_genus"),
                        max_iterations = 10,
                        p_cutoff = 10^-4,
                        score_cutoff = 2,
                        min_peptides = 50,
                        exclusion_method = "genus") {

  # 1. Handle Beads / Peptide Exclusion
  if (!is.null(beads_controls)) {
    # Identify peptides reactive in beads (> 1 FC)
    # We pivot to find any fc > 1 in any bead column
    bead_hits <- beads_controls %>%
      dplyr::select(any_of(c("pep_aa", annotation_cols)), everything()) %>%
      tidyr::pivot_longer(
        cols = -any_of(annotation_cols),
        names_to = "sample",
        values_to = "fc"
      ) %>%
      dplyr::filter(fc > 1) %>%
      dplyr::distinct(pep_aa) %>%
      dplyr::pull(pep_aa)

    bad_peptides <- unique(c(bad_peptides, bead_hits))
  }

  # Filter Data
  clean_fc <- fold_change
  if (!is.null(bad_peptides)) {
    clean_fc <- clean_fc %>%
      dplyr::filter(!pep_aa %in% bad_peptides)
  }

  # 2. Pivot to Long Format and log2 transform fold changes
  # Dynamic pivoting: assume everything NOT an annotation is a sample
  long_data <- clean_fc %>%
    tidyr::pivot_longer(
      cols = -any_of(annotation_cols),
      names_to = "sample_id",
      values_to = "fc"
    ) %>%
    mutate(log2fc = log2(fc)) %>%
    select(-fc)

  # 3. Pre-calculate Group Metrics (Aggregation)
  # Sum scores and count peptides per species/genus/group
  grouped_metrics <- long_data %>%
    dplyr::group_by(taxon_species, taxon_genus, sample_id) %>%
    dplyr::summarise(
      score = sum(log2fc, na.rm = TRUE),
      total_peps = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(score_norm = score / total_peps) %>%
    dplyr::filter(total_peps >= min_peptides)

  # 4. Iterate over Samples
  unique_samples <- unique(grouped_metrics$sample_id)

  # Initialize progressor
  p <- progressr::progressor(steps = length(unique_samples))

  # Using future_map_dfr to loop and automatically bind rows in parallel
  final_results <- furrr::future_map_dfr(unique_samples, function(current_sample) {

    p() # Signal progress

    # Subset data for this sample
    sample_metrics <- grouped_metrics %>% dplyr::filter(sample_id == current_sample)
    sample_peptides <- long_data %>% dplyr::filter(sample_id == current_sample)

    # Run the iterative algorithm
    run_iterative_landscape(
      norm_log = sample_metrics,
      all_peptide_fcs = sample_peptides,
      max_iterations = max_iterations,
      p_cutoff = p_cutoff,
      score_cutoff = score_cutoff,
      exclusion_method = exclusion_method
    )

  }, .options = furrr::furrr_options(seed = TRUE))

  return(final_results)
}

#' Internal Iterative Solver
#'
#' Repeatedly calls calculate_landscape, updating the list of "positive" hits
#' to exclude from the null distribution in the next round.
#'
#' @noRd
run_iterative_landscape <- function(norm_log,
                                    all_peptide_fcs,
                                    max_iterations,
                                    p_cutoff,
                                    score_cutoff,
                                    exclusion_method) {

  iteration <- 0
  # Start with empty positives (nothing excluded)
  # We create a dummy dataframe
  current_positives <- norm_log[0, ]

  # Store history of positives to check for convergence
  history_count <- 0

  scores <- NULL

  while (iteration < max_iterations) {

    # Call the statistical engine
    scores <- calculate_landscape(
      norm_log = norm_log,
      all_peptide_fcs = all_peptide_fcs,
      positives = current_positives,
      exclusion_method = exclusion_method
    )

    # Update criteria for "Positives" (Hits)
    # Criterion 1: Significant p-val AND high score
    hits <- scores %>%
      dplyr::filter(p_val < p_cutoff & ARscore > score_cutoff)

    # Merge unique hits
    new_positives <- hits %>%
      dplyr::distinct(taxon_species, taxon_genus, .keep_all = TRUE)

    # Convergence Check
    # If no NEW positives were added compared to last time, we can stop early
    if (iteration > 0) {
      if (nrow(new_positives) <= history_count) {
        break
      }
    }

    current_positives <- new_positives
    history_count <- nrow(current_positives)
    iteration <- iteration + 1
  }

  # Select final columns to return
  output <- scores %>%
    dplyr::select(
      taxon_species, taxon_genus, sample_id,
      total_peps, score_norm, ARscore, nlog_p
    )

  return(output)
}
