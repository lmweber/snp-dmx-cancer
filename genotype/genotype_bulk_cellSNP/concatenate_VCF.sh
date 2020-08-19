#!/bin/bash

#######################################
# Shell script to concatenate VCF files
#######################################

# This script concatenates VCF files from all 3 samples in the HGSOC dataset, 
# generated with the scripts "genotype_bulk_HGSOC_17667XX.sh".


# qsub -V -cwd -l mem_free=2G,h_vmem=3G,h_fsize=100G concatenate_VCF.sh


# start runtime
start=`date +%s`


# ---------------------
# Concatenate VCF files
# ---------------------

# convert gzipped output files to bgzipped format (required by vcftools)

gunzip -c ../../../genotype/17667X1/cellSNP_bulk/cellSNP.cells.vcf.gz > ../../../genotype/17667X1/cellSNP_bulk/cellSNP.cells-bgz.vcf
bgzip ../../../genotype/17667X1/cellSNP_bulk/cellSNP.cells-bgz.vcf
gunzip -c ../../../genotype/17667X2/cellSNP_bulk/cellSNP.cells.vcf.gz > ../../../genotype/17667X2/cellSNP_bulk/cellSNP.cells-bgz.vcf
bgzip ../../../genotype/17667X2/cellSNP_bulk/cellSNP.cells-bgz.vcf
gunzip -c ../../../genotype/17667X3/cellSNP_bulk/cellSNP.cells.vcf.gz > ../../../genotype/17667X3/cellSNP_bulk/cellSNP.cells-bgz.vcf
bgzip ../../../genotype/17667X3/cellSNP_bulk/cellSNP.cells-bgz.vcf


# concatenate VCF files using vcftools (vcf-concat)

mkdir -p ../../../genotype/cellSNP_bulk_merged

vcf-concat ../../../genotype/17667X1/cellSNP_bulk/cellSNP.cells-bgz.vcf.gz ../../../genotype/17667X2/cellSNP_bulk/cellSNP.cells-bgz.vcf.gz ../../../genotype/17667X3/cellSNP_bulk/cellSNP.cells-bgz.vcf.gz > \
../../../genotype/cellSNP_bulk_merged/cellSNP.cells-merged.vcf


# keep both unzipped and gzipped versions

gzip -c ../../../genotype/cellSNP_bulk_merged/cellSNP.cells-merged.vcf > ../../../genotype/cellSNP_bulk_merged/cellSNP.cells-merged.vcf.gz


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../genotype/runtimes/genotype_bulk_cellSNP
echo runtime: $runtime seconds > ../../../genotype/runtimes/genotype_bulk_cellSNP/runtime_concatenate_VCF.txt

