#!/bin/bash

# ------------------------------------------------
# Shell script to merge and index parsed BAM files
# ------------------------------------------------

# notes:
# - BAM files from Cell Ranger are already position sorted, so do not need to sort
# - there are 3 samples in HGSOC dataset, and 6 samples in lung dataset

# runtime: ~4-6 hours

# qsub -V -cwd -l mem_free=10G,h_vmem=20G,h_fsize=100G merge_and_index_BAM.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: output directory
# $5: sample ID 1
# $6: sample ID 2
# $7: sample ID 3
# $8: sample ID 4
# $9: sample ID 5
# $10: sample ID 6


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


mkdir -p $4/bam_merged

# merge BAM files
samtools merge $4/bam_merged/bam_merged.bam \
$4/$5/outs/possorted_genome_bam_parsed.bam $4/$6/outs/possorted_genome_bam_parsed.bam $4/$7/outs/possorted_genome_bam_parsed.bam \
$4/$8/outs/possorted_genome_bam_parsed.bam $4/$9/outs/possorted_genome_bam_parsed.bam $4/${10}/outs/possorted_genome_bam_parsed.bam

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

