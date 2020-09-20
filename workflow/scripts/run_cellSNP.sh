#!/bin/bash

# ---------------------------
# Shell script to run cellSNP
# ---------------------------

# run cellSNP to genotype cells

# notes:
# - running cellSNP in mode 1
# - using .vcf file from best-performing option for genotyping step (matched bulk 
# RNA-seq samples using bcftools)
# - requires merged BAM file and merged cell barcodes file (each containing unique 
# sample IDs), and .vcf file from genotyping step

# for more details:
# - https://vireosnp.readthedocs.io/en/latest/genotype.html
# - https://github.com/single-cell-genetics/cellSNP

# runtime: ~4 hours (with 10 cores)

# qsub -V -cwd -pe local 10 -l mem_free=5G,h_vmem=10G,h_fsize=100G run_cellSNP.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: genotype directory
# $5: output directory


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------

# note: vcf file needs to be uncompressed
# if still in "vcf.bgz" or "vcf.gz" format then uncompress first
# (for .bgz format, can rename to .gz then gunzip)


cellSNP \
-s $5/bam_merged/bam_merged.bam \
-b $5/barcodes_merged/barcodes_merged.tsv \
-O $5/cellSNP \
-R $4/bcftools/bcftools_HGSOC_rehead.vcf \
-p $3 \
--minMAF=0.05


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/cellSNP
echo runtime: $runtime seconds > $1/cellSNP/runtime_cellSNP.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/cellSNP
date > $2/cellSNP/timestamp_cellSNP.txt
# -----------------------------------

