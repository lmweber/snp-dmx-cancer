#!/bin/bash

##################################################################
# Shell script to run cellSNP and Vireo using single-cell VCF file
##################################################################

# This script runs cellSNP and Vireo on the merged BAM file from the HGSOC 
# dataset, using a custom VCF file generated from genotyping single-cell 
# samples.

# The previous script in this pipeline is "remove_dups.R".


# qsub -V -cwd -pe local 20 -l mem_free=2G,h_vmem=3G,h_fsize=100G run_cellSNP_Vireo_singlecell_VCF.sh


# -----------
# Run cellSNP
# -----------

# note: more stable to run cellSNP interactively using qrsh instead of qsub; not sure why

cellSNP \
-s ../../outputs/HGSOC/bam_merged/bam_merged.bam \
-b ../../outputs/HGSOC/barcodes_merged/barcodes_merged.tsv \
-O ../../genotype_singlecell/cellSNP \
-R ../../genotype_singlecell/cellSNP_singlecell_merged/cellSNP.cells-merged-nodups.vcf \
-p 20 \
--minMAF=0.01


# ---------
# Run Vireo
# ---------

# note: parameter for known number of samples (3 for HGSOC dataset, 6 for lung dataset)

vireo \
-c ../../genotype_singlecell/cellSNP \
-N 3 \
-o ../../genotype_singlecell/vireo \
--randSeed=123

