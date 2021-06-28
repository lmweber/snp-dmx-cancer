###########################################################################
# Debris simulations: generate awk lookup tables and updated barcodes files
###########################################################################

# Simulate debris by assigning reads from some percentage of cell barcodes to
# other cell barcodes in the merged BAM file.

# This script generates:
# (i) lookup tables for each debris scenario, containing the sets of lysed cell
# barcodes and remaining cell barcodes. These lookup tables are then used in the
# awk command in the shell script for the simulation scenario.
# (ii) updated merged cell barcodes file for each scenario


# qrsh -l mem_free=2G,h_vmem=3G,h_fsize=100G
# module load conda_R/4.0
# Rscript generate_awk_lookup_tables_debris.R


# ---------------------
# Simulation parameters
# ---------------------

# parameters for each simulation scenario

# proportion of lysed cells
prop_debris_sims <- c(0.1, 0.2)

# dataset names
dataset_name_sims <- c("HGSOC", "lung")


# ------------------------------------------
# Function to generate lookup tables for awk
# ------------------------------------------

# function to save lookup tables for each simulation scenario

# arguments:
# prop_lysed: proportion of lysed cells to simulate
# dataset_name: dataset name for output files
# file_barcodes_merged: tsv file containing merged list of cell barcodes from Cell Ranger
f_sim_debris <- function(prop_debris, dataset_name, file_barcodes_merged) {
  
  library(tidyverse)
  library(magrittr)
  
  
  # ----------------------
  # Generate lookup tables
  # ----------------------
  
  # load merged list of cell barcodes
  df_barcodes <- read_table(file_barcodes_merged, col_names = "barcode")
  
  # number of cells
  n_cells <- nrow(df_barcodes)
  # number of cells to lyse
  n_lysed <- round(prop_debris * n_cells)
  
  print(n_cells)
  print(n_lysed)
  
  # select random sets of cells to lyse
  set.seed(123)
  ix_lysed <- sample(seq_len(n_cells), n_lysed)
  ix_remaining <- setdiff(seq_len(n_cells), ix_lysed)
  
  # generate lookup table / list of lysed cell barcodes
  df_lookup_lysed <- tibble(
    index = seq_along(ix_lysed), 
    lysed = df_barcodes$barcode[ix_lysed]
  )
  
  # generate lookup table / list of remaining cell barcodes
  df_lookup_remaining <- tibble(
    index = seq_along(ix_remaining), 
    remaining = df_barcodes$barcode[ix_remaining]
  )
  
  # save lookup tables
  fn_out_lysed <- file.path(
    paste0("../../../supplementary_debris/scenarios/", dataset_name, "/30pc"), 
    paste0("debris_lysed_", dataset_name, "_doublets30pc_debris", prop_debris * 100, "pc.tsv")
  )
  fn_out_remaining <- file.path(
    paste0("../../../supplementary_debris/scenarios/", dataset_name, "/30pc"), 
    paste0("debris_remaining_", dataset_name, "_doublets30pc_debris", prop_debris * 100, "pc.tsv")
  )
  write_tsv(df_lookup_lysed, fn_out_lysed)
  # write_tsv(df_lookup_remaining, fn_out_remaining)
  
  
  # ------------------------------
  # Generate updated barcodes file
  # ------------------------------
  
  # save updated barcodes file containing remaining cells
  # fn_out_barcodes <- file.path(
  #   paste0("../../../supplementary_debris/scenarios/", dataset_name, "/30pc"), 
  #   paste0("barcodes_merged_", dataset_name, "_doublets30pc_debris", prop_debris * 100, "pc.tsv")
  # )
  write_tsv(df_lookup_remaining, fn_out_remaining, col_names = FALSE)
  
}


# -----------------------------------------------------------------
# Save lookup tables and barcodes file for each simulation scenario
# -----------------------------------------------------------------

# run function to save lookup tables and updated barcodes file for each
# simulation scenario

for (prop_debris in prop_debris_sims) {
  for (dataset_name in dataset_name_sims) {
    runtime <- system.time({
      file_barcodes_merged <- file.path("../../../benchmarking/scenarios", dataset_name, "30pc", 
                                        paste0("barcodes_merged_", dataset_name, "_30pc.tsv"))
      f_sim_debris(prop_debris, dataset_name, file_barcodes_merged)
    })
  }
}

