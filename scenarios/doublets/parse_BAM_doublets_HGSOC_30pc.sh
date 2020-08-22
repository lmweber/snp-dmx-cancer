#!/bin/bash

#########################################
# Shell script to run doublets simulation
#########################################

# This script runs part of the doublets simulation by:
# (i) converting BAM to SAM
# (ii) parsing through the merged SAM file to replace and combine some 
# percentage of cell barcodes
# (iii) converting SAM back to BAM

# The modifed BAM file can then be used by cellSNP and Vireo (or alternative 
# tools) in the following scripts.

# Notes:
# - lookup table to use in awk command is saved in a .tsv file generated with 
# the script "generate_awk_lookup_tables_doublets.R"
# - for more details on how to use awk for multiple replacements see: 
# https://stackoverflow.com/questions/14234907/replacing-values-in-large-table-using-conversion-table


# qsub -V -cwd -l mem_free=2G,h_vmem=3G,h_fsize=100G parse_BAM_doublets_HGSOC.sh


# ----------------------------------------------------------------
# Generate files containing awk lookup tables and updated barcodes
# ----------------------------------------------------------------

# run R script to generate .tsv files containing awk lookup tables and updated 
# lists of cell barcodes for all simulation scenarios (if not already done)

# qrsh -l mem_free=2G,h_vmem=3G,h_fsize=100G
# module load conda_R/4.0
# Rscript generate_awk_lookup_tables_doublets.R


# -----------------------------
# Parse through merged BAM file
# -----------------------------

# parse through merged BAM file to combine some percentage of cell barcodes

# for one simulation scenario (dataset, percentage of doublets)


# note hyphen for argument order
samtools view -h ../../../scenarios/outputs/HGSOC/bam_merged/bam_merged.bam | \
awk \
'NR==1 { next } FNR==NR { a[$1]=$2; next } (i=gensub(/.*CB\:Z\:([A-Za-z]+\-[A-Za-z0-9]+).*/, "\\1", 1, $0)) in a { gsub(i, a[i]) }1' \
../../../scenarios/doublets/HGSOC/30pc/lookup_table_doublets_HGSOC_30pc.tsv - | \
samtools view -bo ../../../scenarios/doublets/HGSOC/30pc/bam_merged_doublets_HGSOC_30pc.bam


# ---------
# Index BAM
# ---------

samtools index ../../../scenarios/doublets/HGSOC/30pc/bam_merged_doublets_HGSOC_30pc.bam

