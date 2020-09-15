###########################################
# Summary table: number of cells per sample
###########################################


library(tidyverse)
library(magrittr)


# HGSOC dataset

# load ground truth barcodes
fn_barcodes_merged <- "../../benchmarking/outputs/HGSOC/barcodes_merged/barcodes_merged.tsv"
barcodes_merged <- read_table(fn_barcodes_merged, col_names = "barcode_id")

df_truth <- barcodes_merged %>% 
  mutate(sample_id = gsub("^[A-Za-z]+-", "", barcode_id))

table(df_truth$sample_id)


# lung dataset

# load ground truth barcodes
fn_barcodes_merged <- "../../benchmarking/outputs/lung/barcodes_merged/barcodes_merged.tsv"
barcodes_merged <- read_table(fn_barcodes_merged, col_names = "barcode_id")

df_truth <- barcodes_merged %>% 
  mutate(sample_id = gsub("^[A-Za-z]+-", "", barcode_id))

table(df_truth$sample_id)

