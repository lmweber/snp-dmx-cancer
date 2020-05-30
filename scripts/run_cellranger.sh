#!/bin/bash

# script to run Cell Ranger to create BAM files for scRNA-seq samples

# see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct

# notes:
# - maximum file size and memory usage in cluster job submission needs to be large enough, otherwise Cell Ranger fails
# - option '--nosecondary' disables secondary analysis (e.g. dimension reduction) for faster runtime

# runtime: up to 12 hours (using 10 cores)


# arguments:
# $1: sample ID
# $2: directory containing FASTQ files
# $3: directory containing transcriptome reference
# $4: directory to run Cell Ranger (determines output directory)
# $5: directory for timestamp files (for Snakemake)


cwd=$(pwd)
mkdir -p $4
cd $4

cellranger count --id=$1 \
--fastqs=$2 \
--sample=$1 \
--transcriptome=$3 \
--nosecondary \
--jobmode=local \
--localcores=10 \
--localmem=50


# save timestamp file (for Snakemake)
cd $cwd
mkdir -p $5
date > $5/timestamp_cellranger_$1.txt


