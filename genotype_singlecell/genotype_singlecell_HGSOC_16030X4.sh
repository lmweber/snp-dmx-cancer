#!/bin/bash

######################################################
# Shell script to genotype single-cell RNA-seq samples
######################################################

# This script runs cellSNP (mode 2) to generate a custom VCF file by genotyping 
# the single-cell RNA-seq samples (HGSOC dataset). This VCF file can then be 
# provided to cellSNP/Vireo to improve demultiplexing performance.


# qsub -V -cwd -pe local 20 -l mem_free=10G,h_vmem=11G,h_fsize=100G genotype_singlecell.sh


# --------------------
# Run cellSNP (mode 2)
# --------------------

# note: more stable to run interactively with "qrsh" instead of "qsub" on cluster

# gunzip barcodes file
gunzip -c ../../outputs/HGSOC/16030X4_HJTWLDMXX/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ../../genotype_singlecell/16030X4/barcodes.tsv

cellSNP \
-s ../../outputs/HGSOC/16030X4_HJTWLDMXX/outs/possorted_genome_bam.bam \
-b ../../genotype_singlecell/16030X4/barcodes.tsv \
-O ../../genotype_singlecell/16030X4/cellSNP_singlecell \
-p 20 \
--minMAF=0.05 \
--minCOUNT=50 \
--UMItag=None

