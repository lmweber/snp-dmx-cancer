#!/bin/bash

# -------------------------------
# Shell script to run Cell Ranger
# -------------------------------

# for more details:
# - https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct

# notes:
# - maximum file size and memory usage in cluster job submission need to be large enough, otherwise Cell Ranger fails
# - option '--nosecondary' disables secondary analysis (e.g. dimension reduction) for faster runtime

# runtime: ~2-6 hours (with 10 cores)

# qsub -V -cwd -pe local 10 -l mem_free=5G,h_vmem=10G,h_fsize=100G run_cellranger.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: sample ID
# $5: directory containing transcriptome reference
# $6: directory containing FASTQ files for this sample
# $7: output directory


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


# set working directory for output path (Cell Ranger does not have any other 
# option to specify output directory)
cwd=$(pwd)
mkdir -p $7
cd $7


# notes:
# - hard-coding parameters for number of cores, memory, and virtual memory; 
# since other values tend to give unexpected errors on our cluster
cellranger count --id=$4 \
--description=$4 \
--transcriptome=$5 \
--fastqs=$6 \
--sample=$4 \
--nosecondary \
--jobmode=local \
--localcores=10 \
--localmem=50 \
--localvmem=100


# restore working directory
cd $cwd


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/cellranger
echo runtime: $runtime seconds > $1/cellranger/runtime_cellranger_$4.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/cellranger
date > $2/cellranger/timestamp_cellranger_$4.txt
# -----------------------------------

