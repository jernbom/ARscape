# This script is intended to be run on an HPC

# -------------------------------------------------------------------------
# SETUP & LIBRARIES
# -------------------------------------------------------------------------

library(tidyverse)
library(future)
library(ARscape)

options(future.globals.maxSize = 4 * 1024^3) # 4 * 1024^3 = 4 GB

# -------------------------------------------------------------------------
# ARGS
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
input_fc <- args[1]
input_hfc <- args[2]
output_arscores <- args[3]

# -------------------------------------------------------------------------
# LOAD DATA
# -------------------------------------------------------------------------

message("Input paths:")
message(input_fc)
message(input_hfc)

fc <- readr::read_tsv(input_fc, show_col_types = FALSE)
hfc <- readr::read_tsv(input_hfc, show_col_types = FALSE)
mock_ips <- fc %>% colnames() %>% {.[grepl("bead",.,ignore.case = TRUE)]}

# -------------------------------------------------------------------------
# CONFIGURE PARALLEL PLAN
# -------------------------------------------------------------------------

# 1. Detect the number of allocated cores from the environment.
n_cores <- parallelly::availableCores(omit = 1) # Omit 1 to maintain responsiveness when debugging
message(sprintf("System has %d available cores (omitting 1 for stability).", n_cores))

# 2. Set the 'future' plan using all available cores
message(paste("Setting up cluster:"))
cl <- parallel::makeForkCluster(n_cores)
message(paste("A cluster has been set up:"))
print(cl)
plan(cluster, workers = cl)
message("Plan set using above cluster.")

# -------------------------------------------------------------------------
# EXECUTION
# -------------------------------------------------------------------------

t0 <- Sys.time()
message(paste("Started at:", t0))

output <-
  run_arscape(
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
message(paste("Finished at:", t1))
print(t1 - t0)

# -------------------------------------------------------------------------
# SAVE
# -------------------------------------------------------------------------

readr::write_rds(output, output_arscores)
message(paste("Results saved to:", output_arscores))

