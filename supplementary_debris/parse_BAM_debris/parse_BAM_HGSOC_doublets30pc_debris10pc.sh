#!/bin/bash

#######################################
# Shell script to run debris simulation
#######################################

# This script runs part of the debris simulation by:
# (i) converting BAM to SAM
# (ii) parsing through the merged SAM file to lyse some percentage of cell 
# barcodes and assign these reads to other cell barcodes
# (iii) converting SAM back to BAM

# The modifed BAM file can then be used by cellSNP and Vireo (or alternative 
# tools) in the following scripts.

# Notes:
# - lookup tables used in awk command are saved in .tsv files generated with 
# the scripts "generate_awk_lookup_tables_X_debris.R"


# qsub -cwd -l mem_free=2G,h_vmem=3G,h_fsize=100G parse_BAM_debris.sh


# start runtime
start=`date +%s`


# -----------------------------
# Parse through merged BAM file
# -----------------------------

# parse through merged BAM file to lyse some percentage of cell barcodes and 
# assign these reads to other cell barcodes

# for one simulation scenario (dataset, percent doublets, percent debris)


# note hyphen for argument order
samtools view -h ../../../benchmarking/scenarios/HGSOC/30pc/bam_merged_doublets_HGSOC_30pc.bam | \
awk \
'function assign() { cmd = "shuf -n 1 ../../../supplementary_debris/scenarios/HGSOC/30pc/debris_remaining_HGSOC_doublets30pc_debris10pc.tsv"; cmd | getline assigned; close(cmd); return assigned } 
NR==1 { next } FNR==NR { a[$1]=$1; next } (i=gensub(/.*CB\:Z\:([A-Za-z]+\-[A-Za-z0-9]+).*/, "\\1", 1, $0)) in a { gsub(i, assign()) }1' \
../../../supplementary_debris/scenarios/HGSOC/30pc/debris_lysed_HGSOC_doublets30pc_debris10pc.tsv - | \
samtools view -bo ../../../supplementary_debris/scenarios/HGSOC/30pc/bam_merged_HGSOC_doublets30pc_debris10pc.bam


# ---------
# Index BAM
# ---------

samtools index ../../../supplementary_debris/scenarios/HGSOC/30pc/bam_merged_HGSOC_doublets30pc_debris10pc.bam


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../supplementary_debris/scenarios/HGSOC/30pc
echo runtime: $runtime seconds > ../../../supplementary_debris/scenarios/HGSOC/30pc/runtime_parse_BAM_HGSOC_doublets30pc_debris10pc.txt

