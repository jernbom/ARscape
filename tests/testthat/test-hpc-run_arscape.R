test_that("Large scale SLURM job runs correctly", {
  skip_on_cran()
  skip_if_not(Sys.getenv("RUN_HPC_TESTS") == "true")

  # 1. Locate the payload script
  runner_path  <- normalizePath(test_path("testdata", "slurm_runner.sh"))
  payload_path <- normalizePath(test_path("testdata", "slurm_payload.R"))

  # Check if they exist (sanity check)
  if (!file.exists(runner_path)) {
    stop("Could not find slurm_runner.sh at: ", runner_path)
  }
  if (!file.exists(payload_path)) {
    stop("Could not find slurm_payload.R at: ", payload_path)
  }

  # 2. Setup IO paths
  init_time <- format(Sys.time(),"%Y%m%d_%H%M%S")
  job_time <- 120 # minutes
  job_name <- "testthat_hpc_runner"
  log_filename_pattern <- paste0("testthat_%j_", init_time, ".log")
  log_path_pattern     <- file.path(Sys.getenv("ARSCAPE_SANDBOX_PATH"), "log", log_filename_pattern)

  input_fc <- normalizePath(file.path(Sys.getenv("ARSCAPE_TEST_DATA_PATH"), "fiftysixmer_VirscanLar_000_FoldChange_annotated.tsv"))
  input_hfc <- normalizePath(file.path(Sys.getenv("ARSCAPE_TEST_DATA_PATH"), "fiftysixmer_VirscanLar_000_Hits_foldchange_annotated.tsv"))
  output_arscores <- file.path(Sys.getenv("ARSCAPE_SANDBOX_PATH"), "results", paste0("test_hpc_run_arscape_", init_time, ".rds"))

  # 3. Construct the sbatch command
  cmd <- paste(
    "/cm/shared/apps/slurm/current/bin/sbatch",
    "--job-name", job_name,
    "--output", log_path_pattern, # Pass the pattern with %j
    "--time", job_time,
    "--mail-user", Sys.getenv("USER_EMAIL"),

    # The Runner Script
    runner_path,

    # Arguments passed to the Runner Script
    payload_path,
    input_fc,
    input_hfc,
    output_arscores
  )

  # 4. Submit and Capture Job ID
  # We use intern=TRUE to capture the output: "Submitted batch job [NNN]"
  submission_output <- system(cmd, intern = TRUE)

  expect_true(grepl("Submitted batch job", submission_output))

  # Extract the digits from the output string
  job_id <- stringr::str_extract(submission_output, "\\d+")

  if (is.na(job_id)) {
    stop("Failed to parse Job ID. Submission output was: ", submission_output)
  }

  message(paste("SLURM Job ID detected:", job_id))

  # 5. Polling
  job_done <- FALSE
  start_time <- Sys.time()

  # Reconstruct the log filename
  log_filename <- file.path(Sys.getenv("ARSCAPE_SANDBOX_PATH"), "log", paste0("testthat_", job_id, "_", init_time, ".log"))

  while(!job_done) {
    queue_cmd <- paste("/cm/shared/apps/slurm/current/bin/squeue --job", job_id)
    queue_status <- system(queue_cmd, intern = TRUE)

    if (length(queue_status) < 2) { # Header only = job gone (finished)
      job_done <- TRUE
    } else {
      Sys.sleep(30) # Wait 30 seconds before checking again

      # Timeout failsafe
      if (difftime(Sys.time(), start_time, units="mins") > job_time) {
        stop("Test timed out waiting for SLURM job ", job_id)
      }
    }
  }

  # 6. Assess Results

  # A. Check the File
  expect_true(file.exists(output_arscores))

  # Check if file is larger than 50MB
  actual_size <- if(file.exists(output_arscores)) file.size(output_arscores) else 0
  expect_gt(actual_size, 50000000)

  # B. Check the Log (Progress Bar Verification)
  expect_true(file.exists(log_filename))

  if (file.exists(log_filename)) {
    log_lines <- readLines(log_filename)
    # Did the progress bar reach 100%?
    expect_true(any(grepl("Progress: ──────────────────────────────────────────────────────────────── 100%", log_lines)))
  }
})
