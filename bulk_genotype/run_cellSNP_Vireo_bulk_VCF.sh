#!/bin/bash

###########################################################
# Shell script to run cellSNP and Vireo using bulk VCF file
###########################################################

# This script runs cellSNP and Vireo on the merged BAM file from the HGSOC 
# dataset, using a custom VCF file generated from genotyping matched bulk 
# samples.

# The previous script in this pipeline is "remove_dups.R".


# qsub -V -cwd -pe local 20 -l mem_free=2G,h_vmem=3G,h_fsize=100G run_cellSNP_Vireo_bulk_VCF.sh


# -----------
# Run cellSNP
# -----------

# note: more stable to run cellSNP interactively using qrsh instead of qsub; not sure why

cellSNP \
-s ../../outputs/HGSOC/bam_merged/bam_merged.bam \
-b ../../outputs/HGSOC/barcodes_merged/barcodes_merged.tsv \
-O ../../bulk/cellSNP \
-R ../../bulk/cellSNP_merged/cellSNP.cells-merged-nodups.vcf \
-p 20 \
--minMAF=0.01


# ---------
# Run Vireo
# ---------

# note: parameter for known number of samples (3 for HGSOC dataset, 6 for lung dataset)

vireo \
-c ../../bulk/cellSNP \
-N 3 \
-o ../../bulk/vireo \
--randSeed=123

