##################
# Evaluation plots
##################

# Script to generate plots for performance evaluations of benchmarking scenarios

# HGSOC dataset, no doublets


# module load conda_R/4.0
# Rscript evaluations.R


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


# ----------------------------------
# load benchmarking scenario outputs
# ----------------------------------

# Vireo scenarios
# outputs are saved in "donor_ids.tsv" for each scenario

out_vireo_1000GenomesFilt_cellSNPVireo <- read_tsv("../../benchmarking/scenarios/HGSOC/nodoublets/1000GenomesFilt_cellSNPVireo/vireo/donor_ids.tsv")
out_sub <- out_vireo_1000GenomesFilt_cellSNPVireo[, c("cell", "donor_id", "best_doublet")]

# match rows

# no doublets scenarios: same number of rows
stopifnot(nrow(df_truth) == nrow(out_sub))

df_truth <- right_join(df_truth, out_sub, by = c("barcode_id" = "cell"))

df_truth$truth <- factor(df_truth$sample_id)
df_truth$predicted <- factor(df_truth$donor_id)


# ------------------------------
# calculate precision and recall
# ------------------------------

# match sample IDs

# summary table
tbl_summary <- table(truth = df_truth$truth, predicted = df_truth$predicted)
tbl_summary

# manually match sample IDs (note: this needs to be done manually, since
# predicted sample IDs are in arbitrary order)
levels(df_truth$predicted)[1:3] <- c("X3", "X4", "X2")
levels(df_truth$predicted)

tbl_summary <- table(truth = df_truth$truth, predicted = df_truth$predicted)
tbl_summary

