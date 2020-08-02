#!/bin/bash

#########################################
# Shell script to run doublets simulation
#########################################

# This script runs the doublets simulation by:
# (i) converting BAM to SAM
# (ii) parsing through the merged SAM file to replace and combine some 
# percentage of cell barcodes
# (iii) converting SAM back to BAM

# Note: sed commands are saved in a .tsv file generated with the script 
# "generate_sed_cmds_doublets.R" for each simulation scenario.

# Then continue by running cellSNP and Vireo on the modified BAM file.


# -----------------------------------------------------------
# Generate files containing sed commands and updated barcodes
# -----------------------------------------------------------

# run R script to generate .tsv files containing sed commands and updated lists 
# of cell barcodes for all simulation scenarios (if not already done)

# run on JHPCE cluster
#module load conda_R/4.0

#Rscript generate_sed_cmds_doublets.R


# -----------------------------
# Parse through merged BAM file
# -----------------------------

# run doublets simulation by parsing through merged BAM file to combine some 
# percentage of cell barcodes

# for one simulation scenario (dataset, percentage of doublets)

samtools view -h ../../outputs/lung/BAM_merged/BAM_merged.bam | \
sed -f sed_cmds_doublets_lung_8pc.tsv | \
samtools view -bo ../../outputs/lung/doublets/8pc/BAM_merged_doublets_lung_8pc.bam


# -----------
# Run cellSNP
# -----------

cellSNP \
-s ../../outputs/lung/doublets/8pc/BAM_merged_doublets_lung_8pc.bam \
-b ../../outputs/lung/doublets/8pc/barcodes_merged_lung_8pc.tsv \
-O ../../outputs/lung/doublets/8pc/cellSNP \
-R ../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf \
-p 10 \
--minMAF=0.05


# ---------
# Run Vireo
# ---------

# note parameter for known number of samples

vireo \
-c ../../outputs/lung/doublets/8pc/cellSNP \
-N 6 \
-o ../../outputs/lung/doublets/8pc/vireo \
--randSeed=123

