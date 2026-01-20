test_that("Large scale job runs correctly using multicore.", {
  skip_on_cran()
  # Only run if explicitly requested
  skip_if_not(Sys.getenv("RUN_HPC_TESTS") == "true")

  # 1. Setup & Resource Check
  # We require a compute node (detectCores >= 16) to run this heavy test.
  n_cores <- parallelly::availableCores(omit = 1)

  if ((n_cores + 1) < 16) {
    skip("Skipping HPC test: Not enough cores. Are you on a login node? Please run on a compute node.")
  }

  # Load necessary libraries for the test environment
  message(cat("\n"))
  library(tidyverse)
  library(future)
  library(ARscape)

  # 2. IO Paths
  init_time <- format(Sys.time(),"%Y%m%d_%H%M%S")

  # Inputs
  input_fc_path  <- normalizePath(file.path(Sys.getenv("ARSCAPE_TEST_DATA_PATH"), "fiftysixmer_VirscanLar_000_FoldChange_annotated.tsv"))
  input_hfc_path <- normalizePath(file.path(Sys.getenv("ARSCAPE_TEST_DATA_PATH"), "fiftysixmer_VirscanLar_000_Hits_foldchange_annotated.tsv"))

  # Outputs
  output_arscores <- file.path(Sys.getenv("ARSCAPE_SANDBOX_PATH"), "results", paste0("test_hpc_run_arscape_", init_time, ".rds"))

  # 3. Load & Prep Data
  message("Loading input data...")
  fc  <- readr::read_tsv(input_fc_path, show_col_types = FALSE)
  hfc <- readr::read_tsv(input_hfc_path, show_col_types = FALSE)

  # Extract mock IPs
  mock_ips <- fc %>% colnames() %>% {.[grepl("bead", ., ignore.case = TRUE)]}
  expect_gt(length(mock_ips), 0) # Sanity check

  # 4. Configure Parallel Plan
  # Set max globals size to 4GB
  options(future.globals.maxSize = 4 * 1024^3)

  message(cat("\n"))
  message(sprintf("Setting up cluster with %d workers...", n_cores))
  cl <- parallel::makeForkCluster(n_cores)
  plan(cluster, workers = cl)

  # Ensure we clean up the cluster even if the test fails
  on.exit({
    message("Stopping cluster...")
    parallel::stopCluster(cl)
  }, add = TRUE)

  # 5. Execution
  message("Starting run_arscape...")
  t0 <- Sys.time()
  message(t0)

  output <- run_arscape(
    fold_change = fc,
    hits_fold_change = hfc,
    mock_ips = mock_ips,
    excluded_peptides = NULL,
    annotation_cols = c("u_pep_id", "pep_id", "pos_start", "pos_end", "UniProt_acc", "pep_aa", "taxon_genus", "taxon_species", "gene_symbol", "product"),
    max_iterations = 10,
    p_cutoff = 10^-4,
    score_cutoff = 2,
    min_peptides = 50,
    exclusion_method = "genus",
    progress_bar = TRUE
  )

  t1 <- Sys.time()
  duration <- t1 - t0
  message(paste("Calculation finished in:", round(duration, 2), units(duration)))

  # 6. Save & Assess Results

  # Save the file
  readr::write_rds(output, output_arscores)
  expect_true(file.exists(output_arscores))

  # Check file size
  expect_gt(file.size(output_arscores), 50000000) # > 50MB

  # Check object validity directly
  expect_s3_class(output, "data.frame")
  expect_true("ARscore" %in% names(output))
})
