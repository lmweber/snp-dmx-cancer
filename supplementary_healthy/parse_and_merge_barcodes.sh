#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=10G,h_fsize=100G


# --------------------------------------------------
# Shell script to parse and merge cell barcode files
# --------------------------------------------------

# start runtime
start=`date +%s`

mkdir -p ../../supplementary_healthy/outputs/barcodes_merged


# decompress cell barcode files for each sample
gunzip -c ../../supplementary_healthy/outputs/euts1/euts1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ../../supplementary_healthy/outputs/barcodes_merged/barcodes_euts1.tsv
gunzip -c ../../supplementary_healthy/outputs/nufh3/nufh3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ../../supplementary_healthy/outputs/barcodes_merged/barcodes_nufh3.tsv
gunzip -c ../../supplementary_healthy/outputs/babz3/babz3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ../../supplementary_healthy/outputs/barcodes_merged/barcodes_babz3.tsv
gunzip -c ../../supplementary_healthy/outputs/oaqd2/oaqd2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ../../supplementary_healthy/outputs/barcodes_merged/barcodes_oaqd2.tsv
gunzip -c ../../supplementary_healthy/outputs/ieki3/ieki3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ../../supplementary_healthy/outputs/barcodes_merged/barcodes_ieki3.tsv

# add unique sample IDs to cell barcodes for each sample
sed -i "s|\([A-Z]\+\)\-1|\1\-euts1|g" ../../supplementary_healthy/outputs/barcodes_merged/barcodes_euts1.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-nufh3|g" ../../supplementary_healthy/outputs/barcodes_merged/barcodes_nufh3.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-babz3|g" ../../supplementary_healthy/outputs/barcodes_merged/barcodes_babz3.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-oaqd2|g" ../../supplementary_healthy/outputs/barcodes_merged/barcodes_oaqd2.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-ieki3|g" ../../supplementary_healthy/outputs/barcodes_merged/barcodes_ieki3.tsv

# merge files
cat \
../../supplementary_healthy/outputs/barcodes_merged/barcodes_euts1.tsv \
../../supplementary_healthy/outputs/barcodes_merged/barcodes_nufh3.tsv \
../../supplementary_healthy/outputs/barcodes_merged/barcodes_babz3.tsv \
../../supplementary_healthy/outputs/barcodes_merged/barcodes_oaqd2.tsv \
../../supplementary_healthy/outputs/barcodes_merged/barcodes_ieki3.tsv > \
../../supplementary_healthy/outputs/barcodes_merged/barcodes_merged.tsv


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../supplementary_healthy/runtimes/parse_and_merge_barcodes
echo runtime: $runtime seconds > ../../supplementary_healthy/runtimes/parse_and_merge_barcodes/runtime_parse_and_merge_barcodes.txt

