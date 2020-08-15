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

gunzip -c ../../genotype_bulk/17667X1/cellSNP_bulk/cellSNP.cells.vcf.gz > ../../genotype_bulk/17667X1/cellSNP_bulk/cellSNP.cells-bgz.vcf
bgzip ../../genotype_bulk/17667X1/cellSNP_bulk/cellSNP.cells-bgz.vcf
gunzip -c ../../genotype_bulk/17667X2/cellSNP_bulk/cellSNP.cells.vcf.gz > ../../genotype_bulk/17667X2/cellSNP_bulk/cellSNP.cells-bgz.vcf
bgzip ../../genotype_bulk/17667X2/cellSNP_bulk/cellSNP.cells-bgz.vcf
gunzip -c ../../genotype_bulk/17667X3/cellSNP_bulk/cellSNP.cells.vcf.gz > ../../genotype_bulk/17667X3/cellSNP_bulk/cellSNP.cells-bgz.vcf
bgzip ../../genotype_bulk/17667X3/cellSNP_bulk/cellSNP.cells-bgz.vcf


# concatenate VCF files using vcftools (vcr-concat)

mkdir -p ../../genotype_bulk/cellSNP_bulk_merged

vcf-concat ../../genotype_bulk/17667X1/cellSNP_bulk/cellSNP.cells-bgz.vcf.gz ../../genotype_bulk/17667X2/cellSNP_bulk/cellSNP.cells-bgz.vcf.gz ../../genotype_bulk/17667X3/cellSNP_bulk/cellSNP.cells-bgz.vcf.gz > \
../../genotype_bulk/cellSNP_bulk_merged/cellSNP.cells-merged.vcf


# keep both unzipped and gzipped versions

gzip -c ../../genotype_bulk/cellSNP_bulk_merged/cellSNP.cells-merged.vcf > ../../genotype_bulk/cellSNP_bulk_merged/cellSNP.cells-merged.vcf.gz

