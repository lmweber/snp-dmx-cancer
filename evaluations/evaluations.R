##################
# Evaluation plots
##################

# Script to generate plots for performance evaluations of benchmarking scenarios

# HGSOC dataset, 30% doublets


library(tidyverse)
library(magrittr)


# -----------------
# load ground truth
# -----------------

# load ground truth stored in sample ID suffixes in cell barcodes

fn_barcodes_merged <- "../../benchmarking/outputs/HGSOC/barcodes_merged/barcodes_merged.tsv"
barcodes_merged <- read_table(fn_barcodes_merged, col_names = "barcode_id")

df_truth <- barcodes_merged %>% 
  mutate(sample_id = gsub("^[A-Za-z]+-", "", barcode_id))

head(df_truth)

