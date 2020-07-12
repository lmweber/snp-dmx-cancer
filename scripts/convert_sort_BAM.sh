#!/bin/bash

# ---------------------------------------------------
# Shell script to convert SAM to BAM / sort BAM files
# ---------------------------------------------------

# convert SAM to BAM files to reduce disk usage and sort BAM files

# runtime: ~2 hours

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G convert_sort_BAM.sh

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


# convert SAM to BAM
samtools view -S -b $5/$4/alevin_mappings/$4.sam > $5/$4/alevin_mappings/$4.bam

# sort BAM
samtools sort -o $5/$4/alevin_mappings/$4.bam $5/$4/alevin_mappings/$4.bam


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/convert_sort_BAM
echo runtime: $runtime seconds > $1/convert_sort_BAM/runtime_convert_sort_BAM_$4.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/convert_sort_BAM
date > $2/convert_sort_BAM/timestamp_convert_sort_BAM_$4.txt
# -----------------------------------

