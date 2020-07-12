#!/bin/bash

# --------------------------------------------------
# Shell script to parse and merge cell barcode files
# --------------------------------------------------

# runtime: seconds

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G parse_and_merge_barcodes.sh

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

# note: there are 3 samples in this dataset


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


mkdir -p $4/barcodes_merged

# extract cell barcodes from alevin output files for each sample
awk '{print $1}' $4/$5/alevin_output/alevin/featureDump.txt | tail -n +2 > $4/barcodes_merged/barcodes_$8.txt
awk '{print $1}' $4/$6/alevin_output/alevin/featureDump.txt | tail -n +2 > $4/barcodes_merged/barcodes_$9.txt
awk '{print $1}' $4/$7/alevin_output/alevin/featureDump.txt | tail -n +2 > $4/barcodes_merged/barcodes_${10}.txt

# add unique sample IDs to cell barcodes for each sample
sed -i "s|\([A-Z]\+\)|\1|-$8/g" $4/barcodes_merged/barcodes_$8.txt
sed -i "s|\([A-Z]\+\)|\1|-$9/g" $4/barcodes_merged/barcodes_$9.txt
sed -i "s|\([A-Z]\+\)|\1|-${10}/g" $4/barcodes_merged/barcodes_${10}.txt

# merge files
cat $4/barcodes_merged/barcodes_$8.txt $4/barcodes_merged/barcodes_$9.txt $4/barcodes_merged/barcodes_${10}.txt > $4/barcodes_merged/barcodes_merged.txt


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/HGSOC/parse_and_merge_barcodes
echo runtime: $runtime seconds > $1/HGSOC/parse_and_merge_barcodes/runtime_parse_and_merge_barcodes.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/HGSOC/parse_and_merge_barcodes
date > $2/HGSOC/parse_and_merge_barcodes/timestamp_parse_and_merge_barcodes.txt
# -----------------------------------

