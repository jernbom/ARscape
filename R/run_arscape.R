#' Run ARscape Analysis Pipeline
#'
#' The main wrapper function to process PhIP-Seq fold-change data and generate
#' Aggregate Reactivity Scores (ARscores).
#'
#' @section Parallelization:
#' This function leverages the `furrr` package for parallel processing.
#'
#' To enable parallel execution:
#' \preformatted{
#' n_cores <- parallelly::availableCores()
#' cl <- parallel::makeForkCluster(n_cores)
#' future::plan(cluster, workers = cl)
#' }
#'
#' @param fold_change Data frame containing peptide-level fold changes.
#'   Must contain annotation columns (specified in `annotation_cols`) and sample columns.
#' @param mock_ips Optional vector of column names in `hits_fold_change` corresponding to mock IP controls. If provided along with `hits_fold_change`,
#'   peptides with hits in these controls are excluded.
#'   In addition, a character vector of `excluded_peptides` can be provided to replace and/or supplement this peptide exclusion criterion.
#' @param excluded_peptides Optional character vector of UNIQUE peptide identifiers to exclude manually. These peptide identifiers MUST be of the same type as the first element in `annotation_cols`, typically `u_pep_id`.
#' @param annotation_cols Character vector of column names that are NOT samples.
#'   Defaults to standard PhIP-seq annotations in the Larman Lab: `c("u_pep_id", "pep_id", "pos_start", "pos_end", "UniProt_acc", "pep_aa", "taxon_genus", "taxon_species", "gene_symbol", "product")`. NB! The first vector element MUST refer to the column containing UNIQUE peptide identifiers (typically `u_pep_id`).
#' @param max_iterations Integer; maximum number of iterations for the background definition algorithm.
#' @param p_cutoff Numeric; -log10 p-value cutoff for defining positive hits. Default is 4 (-log10(0.0001)).
#' @param score_cutoff Numeric; ARscore cutoff for defining positive hits. Default is 2.
#' @param min_peptides Integer; minimum number of peptides required per group to calculate a score. Default is 50.
#' @param exclusion_method Character; method for excluding reactive peptides during iteration.
#'   Options: "genus", "species", "group".
#' @param progress_bar Logical; If using parallel processing, should a progress bar be displayed? No bar will be displayed if `future::plan(sequential)` is used, see `furrr::future_map()`.
#'
#' @return A data frame containing the final ARscores and p-values for all samples.
#'
#' @importFrom dplyr select filter group_by summarise mutate left_join anti_join pick everything all_of any_of n distinct bind_rows pull
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map_dfr
#' @importFrom furrr future_pmap_dfr furrr_options
#' @export
run_arscape <- function(fold_change,
                        hits_fold_change = NULL,
                        mock_ips = NULL,
                        excluded_peptides = NULL,
                        annotation_cols = c("u_pep_id", "pep_id", "pos_start", "pos_end", "UniProt_acc", "pep_aa", "taxon_genus", "taxon_species", "gene_symbol", "product"),
                        max_iterations = 10,
                        p_cutoff = 10^-4,
                        score_cutoff = 2,
                        min_peptides = 50,
                        exclusion_method = "genus",
                        progress_bar = FALSE) {

  # 1. Handle Beads / Peptide Exclusion
  if (!is.null(mock_ips) && !is.null(hits_fold_change)) {
    # Identify peptides reactive in mock IPs.
    mock_ip_hits <-
      hits_fold_change %>%
      dplyr::select(tidyselect::all_of(c(annotation_cols, mock_ips))) %>%
      tidyr::pivot_longer(
        cols = -tidyselect::any_of(annotation_cols),
        names_to = "sample_id",
        values_to = "hfc"
      ) %>%
      dplyr::filter(hfc > 1) %>%
      dplyr::pull(annotation_cols[[1]]) %>%
      unique()

    excluded_peptides <- unique(c(excluded_peptides, mock_ip_hits))
  }


  # Filter Data
  clean_fc <- fold_change
  if (!is.null(excluded_peptides)) {
    clean_fc <- clean_fc %>%
      dplyr::filter(!.data[[annotation_cols[[1]]]] %in% excluded_peptides)
  }

  # 2. Pivot to Long Format and log2 transform fold changes
  long_data <- clean_fc %>%
    tidyr::pivot_longer(
      cols = -tidyselect::all_of(annotation_cols),
      names_to = "sample_id",
      values_to = "fc"
    ) %>%
    mutate(log2fc = log2(fc) %>% {if_else(. > 0, ., 0)}) %>% # Log2FC floored to 0.
    select(-fc)

  # 3. Pre-calculate Group Metrics (Aggregation)
  # Sum scores and count peptides per species/genus/group
  grouped_metrics <- long_data %>%
    dplyr::group_by(taxon_species, taxon_genus, sample_id) %>% # Make taxon_species, taxon_genus dynamic if needed later
    dplyr::summarise(
      score = sum(log2fc, na.rm = TRUE),
      total_peps = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(score_norm = score / total_peps) %>%
    dplyr::filter(total_peps >= min_peptides)

  # 4. Prepare for Parallel Execution
  # Identify valid species (those that met the min_peptides criteria)
  valid_species <-
    grouped_metrics %>%
    dplyr::transmute(genus_species = paste0(taxon_genus, "_", taxon_species)) %>%
    dplyr::pull(genus_species) %>%
    unique()

  if (length(valid_species) == 0) {
    warning("No species met the minimum peptide criteria.")
    return(data.frame())
  }

  # Filter long_data to only valid species to ensure alignment
  long_data <-
    long_data %>%
    dplyr::mutate(genus_species = paste0(taxon_genus, "_", taxon_species)) %>%
    dplyr::filter(genus_species %in% valid_species) %>%
    dplyr::select(-genus_species)

  grouped_metrics_list <- split(grouped_metrics, grouped_metrics$sample_id)
  long_data_list <- split(long_data, long_data$sample_id)

  rm(fold_change, hits_fold_change, clean_fc, long_data, grouped_metrics)
  gc()
  print(lobstr::obj_sizes(grouped_metrics_list, long_data_list))

  # 5. Parallel Map
  final_results <- furrr::future_pmap(
    .l = list(
      norm_log = grouped_metrics_list,
      all_peptide_fcs = long_data_list,
      current_sample_id = names(long_data_list)
      ),
    .f = run_iterative_landscape,
    max_iterations = max_iterations,
    p_cutoff = p_cutoff,
    score_cutoff = score_cutoff,
    exclusion_method = exclusion_method,
    # .id = "sample_id",
    .progress = progress_bar,
    .options = furrr::furrr_options(
      seed = 120,
      globals = c("calculate_landscape"),
      packages = c("stats", "fitdistrplus", "dplyr", "limma", "ARscape")
    )
  )

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
                                    current_sample_id,
                                    max_iterations,
                                    p_cutoff,
                                    score_cutoff,
                                    exclusion_method) {

  message(paste("\nProcessing:", current_sample_id))

  iteration <- 0
  current_positives <- norm_log[0, ]

  # Store history of positives to check for convergence
  history_count <- 0

  scores <- NULL

  hits_log_debug <- list()

  while (iteration < max_iterations) {
    message(iteration)

    # Call the statistical engine
    results_debug <- calc_arscore(
      norm_log = norm_log,
      all_peptide_fcs = all_peptide_fcs,
      sample_id = current_sample_id,
      positives = current_positives,
      exclusion_method = exclusion_method
    )
    scores <- results_debug[[1]]

    # Update detected "Positives" (Hits)
    hits <- scores %>%
      dplyr::filter(
        (p_val < p_cutoff & ARscore > score_cutoff) |
          (p_val < 10^-15)
      )

    hits_log_debug[[(iteration+1)]] <- hits

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


  output_debug <- list()
  output_debug[[1]] <- scores
  output_debug[[2]] <- hits_log_debug
  output_debug[[3]] <- results_debug[[2]]
  output_debug[[4]] <- results_debug[[3]]

  return(output_debug)
}
