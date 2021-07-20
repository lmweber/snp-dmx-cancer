#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=3G,h_fsize=100G


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


# start runtime
start=`date +%s`


# -----------------------------
# Parse through merged BAM file
# -----------------------------

# parse through merged BAM file to lyse some percentage of cell barcodes and 
# assign these reads to other cell barcodes

# for one simulation scenario (dataset, percent doublets, percent debris)


# note hyphen for argument order
samtools view -h ../../../benchmarking/scenarios/lung/30pc/bam_merged_doublets_lung_30pc.bam | \
awk \
-v f_remaining="../../../supplementary_debris/scenarios/lung/30pc/debris_remaining_lung_doublets30pc_debris20pc.tsv" \
-v n_remaining="$(wc -l ../../../supplementary_debris/scenarios/lung/30pc/debris_remaining_lung_doublets30pc_debris20pc.tsv | cut -f1 -d' ')" \
'NR==1 { next } 
FNR==NR { lysed[$2]=$2; next } 
FILENAME==f_remaining { remaining[$1]=$2; next } 
(i=gensub(/.*CB\:Z\:([A-Za-z]+\-[A-Za-z0-9]+).*/, "\\1", 1, $0)) in lysed { gsub(i, remaining[int(rand()*n_remaining+1)]) }1' \
../../../supplementary_debris/scenarios/lung/30pc/debris_lysed_lung_doublets30pc_debris20pc.tsv ../../../supplementary_debris/scenarios/lung/30pc/debris_remaining_lung_doublets30pc_debris20pc.tsv - | \
samtools view -bo ../../../supplementary_debris/scenarios/lung/30pc/bam_merged_lung_doublets30pc_debris20pc.bam


# ---------
# Index BAM
# ---------

samtools index ../../../supplementary_debris/scenarios/lung/30pc/bam_merged_lung_doublets30pc_debris20pc.bam


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../supplementary_debris/scenarios/lung/30pc
echo runtime: $runtime seconds > ../../../supplementary_debris/scenarios/lung/30pc/runtime_parse_BAM_lung_doublets30pc_debris20pc.txt

