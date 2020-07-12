#!/bin/bash

# ------------------------------------------------
# Shell script to merge and index parsed BAM files
# ------------------------------------------------

# runtime: ~3 hours

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G merge_and_index_parsed_BAM.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: output directory
# $5: sample ID 1
# $6: sample ID 2
# $7: sample ID 3

# note: there are 3 samples in this dataset


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


# merge BAM files
mkdir -p $4/bam_merged
samtools merge $4/bam_merged/bam_merged.bam $4/$5/alevin_mappings/$5.bam $4/$6/alevin_mappings/$6.bam $4/$7/alevin_mappings/$7.bam

# index merged BAM
samtools index $4/bam_merged/bam_merged.bam


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/merge_and_index_BAM
echo runtime: $runtime seconds > $1/merge_and_index_BAM/runtime_merge_and_index_BAM.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/merge_and_index_BAM
date > $2/merge_and_index_BAM/timestamp_merge_and_index_BAM.txt
# -----------------------------------

