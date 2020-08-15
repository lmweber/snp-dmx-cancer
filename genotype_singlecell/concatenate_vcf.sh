#!/bin/bash

#######################################
# Shell script to concatenate VCF files
#######################################

# This script concatenates VCF files from all 3 samples in the HGSOC dataset, 
# generated with the scripts "genotype_singlecell_HGSOC_16030XX.sh".


# qsub -V -cwd -pe local 20 -l mem_free=2G,h_vmem=3G,h_fsize=100G concatenate_vcf.sh


# ---------------------
# Concatenate VCF files
# ---------------------

# convert gzipped output files to bgzipped format (required by vcftools)

gunzip -c ../../genotype_singlecell/16030X2/cellSNP/cellSNP.cells.vcf.gz > ../../genotype_singlecell/16030X2/cellSNP/cellSNP.cells-bgz.vcf
bgzip ../../genotype_singlecell/16030X2/cellSNP/cellSNP.cells-bgz.vcf
gunzip -c ../../genotype_singlecell/16030X3/cellSNP/cellSNP.cells.vcf.gz > ../../genotype_singlecell/16030X3/cellSNP/cellSNP.cells-bgz.vcf
bgzip ../../genotype_singlecell/16030X3/cellSNP/cellSNP.cells-bgz.vcf
gunzip -c ../../genotype_singlecell/16030X4/cellSNP/cellSNP.cells.vcf.gz > ../../genotype_singlecell/16030X4/cellSNP/cellSNP.cells-bgz.vcf
bgzip ../../genotype_singlecell/16030X4/cellSNP/cellSNP.cells-bgz.vcf


# concatenate VCF files using vcftools (vcr-concat)

mkdir -p ../../genotype_singlecell/cellSNP_singlecell_merged

vcf-concat ../../genotype_singlecell/16030X2/cellSNP/cellSNP.cells-bgz.vcf.gz ../../genotype_singlecell/16030X3/cellSNP/cellSNP.cells-bgz.vcf.gz ../../genotype_singlecell/16030X4/cellSNP/cellSNP.cells-bgz.vcf.gz > \
../../genotype_singlecell/cellSNP_singlecell_merged/cellSNP.cells-merged.vcf


# keep both unzipped and gzipped versions

gzip -c ../../genotype_singlecell/cellSNP_singlecell_merged/cellSNP.cells-merged.vcf > ../../genotype_singlecell/cellSNP_singlecell_merged/cellSNP.cells-merged.vcf.gz

