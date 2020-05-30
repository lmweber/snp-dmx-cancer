#!/bin/bash

# script to run Cell Ranger to create BAM files for scRNA-seq samples

# see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct

# notes:
# - maximum file sizes and memory usage need to be specified correctly when running on Linux cluster
# - using option '--nosecondary' to disable secondary analysis (e.g. dimension reduction) for faster runtime

# runtime: up to 12 hours (using 10 cores)


# arguments:
# $1: sample ID
# $2: directory containing FASTQ files
# $3: directory containing transcriptome reference
# $4: maximum number of cores
# $5: maximum memory (GB)
# $6: original working directory
# $7: directory for timestamp files (for Snakemake)


cellranger count --id=$1 \
--fastqs=$2 \
--sample=$1 \
--transcriptome=$3 \
--nosecondary \
--localcores=$4 \
--localmem=$5


# save timestamp file (for Snakemake)
cd $6
mkdir -p $7
date > $7/timestamp_cellranger_$1.txt


