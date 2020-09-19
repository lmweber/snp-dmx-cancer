#!/bin/bash

# --------------------------------------------------
# Shell script to parse and merge cell barcode files
# --------------------------------------------------

# notes:
# - parse cell barcode files to contain unique sample IDs matching BAM files
# - there are 3 samples in HGSOC dataset, and 6 samples in lung dataset

# runtime: seconds

# qsub -V -cwd -l mem_free=10G,h_vmem=20G,h_fsize=100G parse_and_merge_barcodes.sh

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
# $11: short sample ID 1
# $12: short sample ID 2
# $13: short sample ID 3
# $14: short sample ID 4
# $15: short sample ID 5
# $16: short sample ID 6


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


mkdir -p $4/barcodes_merged

# decompress cell barcode files for each sample
gunzip -c $4/$5/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $4/barcodes_merged/barcodes_${11}.tsv
gunzip -c $4/$6/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $4/barcodes_merged/barcodes_${12}.tsv
gunzip -c $4/$7/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $4/barcodes_merged/barcodes_${13}.tsv
gunzip -c $4/$8/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $4/barcodes_merged/barcodes_${14}.tsv
gunzip -c $4/$9/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $4/barcodes_merged/barcodes_${15}.tsv
gunzip -c $4/${10}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $4/barcodes_merged/barcodes_${16}.tsv

# add unique sample IDs to cell barcodes for each sample
sed -i "s|\([A-Z]\+\)\-1|\1\-${11}|g" $4/barcodes_merged/barcodes_${11}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${12}|g" $4/barcodes_merged/barcodes_${12}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${13}|g" $4/barcodes_merged/barcodes_${13}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${14}|g" $4/barcodes_merged/barcodes_${14}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${15}|g" $4/barcodes_merged/barcodes_${15}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${16}|g" $4/barcodes_merged/barcodes_${16}.tsv

# merge files
cat $4/barcodes_merged/barcodes_${11}.tsv $4/barcodes_merged/barcodes_${12}.tsv $4/barcodes_merged/barcodes_${13}.tsv \
$4/barcodes_merged/barcodes_${14}.tsv $4/barcodes_merged/barcodes_${15}.tsv $4/barcodes_merged/barcodes_${16}.tsv > \
$4/barcodes_merged/barcodes_merged.tsv


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/parse_and_merge_barcodes
echo runtime: $runtime seconds > $1/parse_and_merge_barcodes/runtime_parse_and_merge_barcodes.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/parse_and_merge_barcodes
date > $2/parse_and_merge_barcodes/timestamp_parse_and_merge_barcodes.txt
# -----------------------------------

