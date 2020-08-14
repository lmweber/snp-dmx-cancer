#!/bin/bash

#######################################
# Shell script to concatenate VCF files
#######################################

# This script concatenates VCF files from all 3 samples in the HGSOC dataset, 
# generated with the scripts "genotype_bulk_HGSOC_17667XX.sh".


# qsub -V -cwd -pe local 20 -l mem_free=2G,h_vmem=3G,h_fsize=100G concatenate_vcf.sh


# ---------------------
# Concatenate VCF files
# ---------------------

# convert gzipped output files to bgzipped format (required by vcftools)

gunzip -c ../../bulk/17667X1/cellSNP/cellSNP.cells.vcf.gz > ../../bulk/17667X1/cellSNP/cellSNP.cells-bgz.vcf
bgzip ../../bulk/17667X1/cellSNP/cellSNP.cells-bgz.vcf
gunzip -c ../../bulk/17667X2/cellSNP/cellSNP.cells.vcf.gz > ../../bulk/17667X2/cellSNP/cellSNP.cells-bgz.vcf
bgzip ../../bulk/17667X2/cellSNP/cellSNP.cells-bgz.vcf
gunzip -c ../../bulk/17667X3/cellSNP/cellSNP.cells.vcf.gz > ../../bulk/17667X3/cellSNP/cellSNP.cells-bgz.vcf
bgzip ../../bulk/17667X3/cellSNP/cellSNP.cells-bgz.vcf


# concatenate VCF files using vcftools (vcr-concat)

mkdir -p ../../bulk/cellSNP_bulk_merged

vcf-concat ../../bulk/17667X1/cellSNP/cellSNP.cells-bgz.vcf.gz ../../bulk/17667X2/cellSNP/cellSNP.cells-bgz.vcf.gz ../../bulk/17667X3/cellSNP/cellSNP.cells-bgz.vcf.gz > \
../../bulk/cellSNP_bulk_merged/cellSNP.cells-merged.vcf


# keep both unzipped and gzipped versions

gzip -c ../../bulk/cellSNP_bulk_merged/cellSNP.cells-merged.vcf > ../../bulk/cellSNP_bulk_merged/cellSNP.cells-merged.vcf.gz

