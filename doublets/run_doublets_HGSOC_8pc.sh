#!/bin/bash

#########################################
# Shell script to run doublets simulation
#########################################

# This script runs the doublets simulation by:
# (i) converting BAM to SAM
# (ii) parsing through the merged SAM file to replace and combine some 
# percentage of cell barcodes
# (iii) converting SAM back to BAM
# (iv) then continue by running cellSNP and Vireo on the modified BAM file

# Notes:
# - lookup table to use in awk command is saved in a .tsv file generated with 
# the script "generate_awk_lookup_tables_doublets.R"
# - for more details on how to use awk for multiple replacements see: 
# https://stackoverflow.com/questions/14234907/replacing-values-in-large-table-using-conversion-table


# qsub -V -cwd -pe local 10 -l mem_free=5G,h_vmem=10G,h_fsize=100G run_doublets.sh


# ----------------------------------------------------------------
# Generate files containing awk lookup tables and updated barcodes
# ----------------------------------------------------------------

# run R script to generate .tsv files containing awk lookup tables and updated 
# lists of cell barcodes for all simulation scenarios (if not already done)

# run on JHPCE cluster
#module load conda_R/4.0

#Rscript generate_awk_lookup_tables_doublets.R


# -----------------------------
# Parse through merged BAM file
# -----------------------------

# run doublets simulation by parsing through merged BAM file to combine some 
# percentage of cell barcodes

# for one simulation scenario (dataset, percentage of doublets)


# note hyphen for argument order
samtools view -h ../../outputs/HGSOC/bam_merged/bam_merged.bam | \
awk \
'NR==1 { next } FNR==NR { a[$1]=$2; next } (i=gensub(/.*CB\:Z\:([A-Za-z]+\-[A-Za-z0-9]+).*/, "\\1", 1, $0)) in a { gsub(i, a[i]) }1' \
../../outputs/HGSOC/doublets/8pc/lookup_table_doublets_HGSOC_8pc.tsv - | \
samtools view -bo ../../outputs/HGSOC/doublets/8pc/bam_merged_doublets_HGSOC_8pc.bam


# -----------
# Run cellSNP
# -----------

cellSNP \
-s ../../outputs/HGSOC/doublets/8pc/bam_merged_doublets_HGSOC_8pc.bam \
-b ../../outputs/HGSOC/doublets/8pc/barcodes_merged_HGSOC_8pc.tsv \
-O ../../outputs/HGSOC/doublets/8pc/cellSNP \
-R ../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf \
-p 10 \
--minMAF=0.05


# ---------
# Run Vireo
# ---------

# note parameter for known number of samples (3 for HGSOC dataset, 6 for lung dataset)

vireo \
-c ../../outputs/HGSOC/doublets/8pc/cellSNP \
-N 3 \
-o ../../outputs/HGSOC/doublets/8pc/vireo \
--randSeed=123

