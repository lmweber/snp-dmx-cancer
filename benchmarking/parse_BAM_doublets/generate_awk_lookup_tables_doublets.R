#############################################################################
# Doublets simulations: generate awk lookup tables and updated barcodes files
#############################################################################

# Simulate doublets by combining some percentage of cell barcodes in the merged
# BAM file.

# This script generates:
# (i) a lookup table for each simulation scenario, containing the set of cell
# barcodes to replace and combine. This lookup table is then used in the awk
# command in the shell script for the simulation scenario.
# (ii) updated merged cell barcodes file for each scenario


# qrsh -l mem_free=2G,h_vmem=3G,h_fsize=100G
# module load conda_R/4.0
# Rscript generate_awk_lookup_tables_doublets.R


# ---------------------
# Simulation parameters
# ---------------------

# parameters for each simulation scenario

# proportion of doublets to simulate (i.e. proportion of final number of cells)
prop_doublets_sims <- c(0.2, 0.3)

# dataset names
dataset_name_sims <- c("HGSOC", "lung")


# ------------------------------------------
# Function to generate lookup tables for awk
# ------------------------------------------

# function to save lookup table for each simulation scenario

# arguments:
# prop_doublets: proportion of doublets to simulate
# dataset_name: dataset name for output files
# file_barcodes_merged: tsv file containing merged list of cell barcodes from Cell Ranger
f_sim_doublets <- function(prop_doublets, dataset_name, file_barcodes_merged) {
  
  library(tidyverse)
  library(magrittr)
  
  # corrected proportion of doublets to simulate, i.e. number of cells to
  # replace after taking into account reduction in final number of cells
  prop_doublets_corrected <- prop_doublets / (1 + prop_doublets)
  print(prop_doublets_corrected)
  
  
  # ---------------------
  # Generate lookup table
  # ---------------------
  
  # load merged list of cell barcodes
  df_barcodes <- read_table(file_barcodes_merged, col_names = "barcode")
  
  # number of cells
  n_cells <- nrow(df_barcodes)
  # number of cells to combine with other cells as doublets
  n_doublets <- round(prop_doublets_corrected * n_cells)
  
  print(n_cells)
  print(n_doublets)
  
  # select random sets of (i) cells to replace and (ii) replacements
  set.seed(123)
  ix_original <- sample(seq_len(n_cells), n_doublets)
  ix_replacement <- sample(setdiff(seq_len(n_cells), ix_original), n_doublets)
  
  # generate lookup table / replacement table for cell barcodes
  df_lookup <- tibble(
    original = df_barcodes$barcode[ix_original], 
    replacement = df_barcodes$barcode[ix_replacement]
  )
  
  # save lookup table
  fn_out <- file.path(
    paste0("../../../benchmarking/scenarios/", dataset_name, "/", prop_doublets * 100, "pc"), 
    paste0("lookup_table_doublets_", dataset_name, "_", prop_doublets * 100, "pc.tsv")
  )
  write_tsv(df_lookup, fn_out)
  
  
  # ------------------------------
  # Generate updated barcodes file
  # ------------------------------
  
  # generate updated list of cell barcodes by removing replaced cells
  barcodes_merged_new <- df_barcodes$barcode[!(df_barcodes$barcode %in% df_lookup$original)]
  stopifnot(all(barcodes_merged_new %in% df_barcodes$barcode))
  stopifnot(length(barcodes_merged_new) == n_cells - n_doublets)
  
  # save barcodes file
  fn_out <- file.path(
    paste0("../../../benchmarking/scenarios/", dataset_name, "/", prop_doublets * 100, "pc"), 
    paste0("barcodes_merged_", dataset_name, "_", prop_doublets * 100, "pc.tsv")
  )
  write_tsv(tibble(barcodes_merged_new), fn_out, col_names = FALSE)
  
}


# ----------------------------------------------------------------
# Save lookup table and barcodes file for each simulation scenario
# ----------------------------------------------------------------

# run function to save lookup table and updated barcodes file for each
# simulation scenario

for (prop_doublets in prop_doublets_sims) {
  for (dataset_name in dataset_name_sims) {
    runtime <- system.time({
      file_barcodes_merged <- file.path("../../../benchmarking/outputs", dataset_name, "barcodes_merged/barcodes_merged.tsv")
      f_sim_doublets(prop_doublets, dataset_name, file_barcodes_merged)
    })
    
    # save runtime
    fn_runtime <- file.path(
      paste0("../../../benchmarking/scenarios/", dataset_name, "/", prop_doublets * 100, "pc"), 
      paste0("runtime_lookup_table_doublets_", dataset_name, "_", prop_doublets * 100, "pc.txt")
    )
    sink(fn_runtime)
    cat(paste0("runtime: ", round(runtime[["elapsed"]]), " seconds"))
    sink()
  }
}

