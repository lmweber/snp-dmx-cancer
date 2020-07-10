#!/bin/bash

# ----------------------------------------
# Shell script to convert SAM to BAM files
# ----------------------------------------

# convert SAM to BAM files to reduce disk usage

# runtime: ~1 hour

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G convert_BAM.sh

# arguments:
# $1: sample ID
# $2: directory for runtimes
# $3: directory for timestamp files
# $4: number of threads
# $5: SAM/BAM output directory from salmon alevin


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


samtools view -S -b $5/$1.sam > $5/$1.bam


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $2/HGSOC/convert_BAM
echo runtime: $runtime seconds > $2/HGSOC/convert_BAM/runtime_convert_BAM_$1.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $3/HGSOC/convert_BAM
date > $3/HGSOC/convert_BAM/timestamp_convert_BAM_$1.txt
# -----------------------------------

