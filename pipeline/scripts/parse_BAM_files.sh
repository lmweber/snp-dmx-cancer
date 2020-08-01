#!/bin/bash

# ------------------------------------------------------------------
# Shell script to parse BAM files to add sample IDs to cell barcodes
# ------------------------------------------------------------------

# notes:
# - convert BAM to SAM, parse SAM, then convert back to BAM to save space
# - syntax to search and replace: sed -i "s/regexp/replacement/g"
# - using double quotes allows variables inside sed expression
# - using alternative separator | allows slashes in variable inside sed expression
# - regular expression matches cell barcode with format "CB:Z:AACTTTCAGCGCTCCA-1"
# - we replace the sample suffix "-1" with a unique sample ID, e.g. "-X2"

# runtime: ~1-2 hours

# qsub -V -cwd -l mem_free=10G,h_vmem=20G,h_fsize=100G parse_BAM_files.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: sample ID
# $5: short sample ID
# $6: output directory


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


samtools view -h $6/$4/outs/possorted_genome_bam.bam | \
sed "s|\(CB\:Z\:[A-Z]\+\)\-1|\1\-$5|g" | \
samtools view -bo $6/$4/outs/possorted_genome_bam_parsed.bam


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/parse_BAM_files
echo runtime: $runtime seconds > $1/parse_BAM_files/runtime_parse_BAM_files_$4.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/parse_BAM_files
date > $2/parse_BAM_files/timestamp_parse_BAM_files_$4.txt
# -----------------------------------

