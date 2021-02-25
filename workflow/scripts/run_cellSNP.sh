#!/bin/bash

# ---------------------------
# Shell script to run cellSNP
# ---------------------------

# run cellSNP to genotype cells

# notes:
# - running cellSNP in mode 1
# - using .vcf file from best-performing option for genotyping step (matched bulk 
# RNA-seq samples using bcftools)
# - requires merged BAM file and merged cell barcodes file from previous steps 
# (doublets simulation scenario), and .vcf file from genotyping step

# for more details:
# - https://vireosnp.readthedocs.io/en/latest/genotype.html
# - https://github.com/single-cell-genetics/cellsnp-lite

# runtime: 

# qsub -V -cwd -pe local 10 -l mem_free=5G,h_vmem=10G,h_fsize=100G run_cellSNP.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: output directory
# $5: genotype directory
# $6: dataset name for simulation scenario
# $7: percentage of doublets for simulation scenario (formatted as e.g. "20pc")


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------

# note: vcf file needs to be uncompressed
# if still in "vcf.bgz" or "vcf.gz" format then uncompress first
# (for .bgz format, can rename to .gz then gunzip)

cellsnp-lite \
-s $4/$6/doublets_sims/$7/bam_merged_doublets_$6_$7.bam \
-b $4/$6/doublets_sims/$7/barcodes_merged_$6_$7.tsv \
-O $4/$6/doublets_sims/$7/cellSNP \
-R $5/bcftools/bcftools_HGSOC_rehead.vcf \
-p $3 \
--minMAF=0.1 \
--minCOUNT=20 \
--gzip


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/$6/doublets_sims/$7/cellSNP
echo runtime: $runtime seconds > $1/$6/doublets_sims/$7/cellSNP/runtime_cellSNP_$6_$7.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/$6/doublets_sims/$7/cellSNP
date > $2/$6/doublets_sims/$7/cellSNP/timestamp_cellSNP_$6_$7.txt
# -----------------------------------

