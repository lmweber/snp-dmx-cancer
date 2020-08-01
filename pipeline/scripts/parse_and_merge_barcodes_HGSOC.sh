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
# $8: short sample ID 1
# $9: short sample ID 2
# $10: short sample ID 3


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


mkdir -p $4/barcodes_merged

# decompress cell barcode files for each sample
gunzip -c $4/$5/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $4/barcodes_merged/barcodes_$8.tsv
gunzip -c $4/$6/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $4/barcodes_merged/barcodes_$9.tsv
gunzip -c $4/$7/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $4/barcodes_merged/barcodes_${10}.tsv

# add unique sample IDs to cell barcodes for each sample
sed -i "s|\([A-Z]\+\)\-1|\1\-$8|g" $4/barcodes_merged/barcodes_$8.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-$9|g" $4/barcodes_merged/barcodes_$9.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${10}|g" $4/barcodes_merged/barcodes_${10}.tsv

# merge files
cat $4/barcodes_merged/barcodes_$8.tsv $4/barcodes_merged/barcodes_$9.tsv $4/barcodes_merged/barcodes_${10}.tsv > \
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

