#!/bin/bash

#############################################################
# Shell script to genotype bulk RNA-seq samples using cellSNP
#############################################################

# This script runs cellSNP (mode 2) to genotype and generate VCF files for the 
# bulk RNA-seq samples in our HGSOC dataset.


# note: requires BAM files from previous scripts "align_bulk_HGSOC_17667XX.sh"

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=11G,h_fsize=100G genotype_bulk.sh


# --------------------
# Run cellSNP (mode 2)
# --------------------

# note: can be more stable to run interactively with "qrsh" instead of "qsub" on cluster

cellSNP \
-s ../../../genotype/17667X1/STAR/Aligned.sortedByCoord.out.bam \
-O ../../../genotype/17667X1/cellSNP_bulk \
-p 10 \
--minMAF=0.01 \
--minCOUNT=50 \
--UMItag=None
