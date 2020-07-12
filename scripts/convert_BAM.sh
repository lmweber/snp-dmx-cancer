#!/bin/bash

# ----------------------------------------
# Shell script to convert SAM to BAM files
# ----------------------------------------

# convert SAM to BAM files to reduce disk usage

# runtime: ~1 hour

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G convert_BAM.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: sample ID
# $5: output directory


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


samtools view -S -b $5/$4/alevin_mappings/$4.sam > $5/$4/alevin_mappings/$4.bam


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/convert_BAM
echo runtime: $runtime seconds > $1/convert_BAM/runtime_convert_BAM_$4.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/convert_BAM
date > $2/convert_BAM/timestamp_convert_BAM_$4.txt
# -----------------------------------

