# Parse and merge cell barcodes files (for cellSNP)


# parse cell barcodes from Cell Ranger to contain unique sample IDs (to match 
# merged BAM file), then merge into a single file


#bash parse_and_merge_barcodes.sh


mkdir barcodes
cd barcodes
cp ../16030X2_HJVMLDMXX/outs/filtered_feature_bc_matrix/barcodes.tsv.gz barcodes-X2.tsv.gz
cp ../16030X3_HJTWLDMXX/outs/filtered_feature_bc_matrix/barcodes.tsv.gz barcodes-X3.tsv.gz
cp ../16030X4_HJTWLDMXX/outs/filtered_feature_bc_matrix/barcodes.tsv.gz barcodes-X4.tsv.gz
gunzip *
sed -i "s/\([A-Z]\+\)\-1/\1\-X2/g" barcodes-X2.tsv
sed -i "s/\([A-Z]\+\)\-1/\1\-X3/g" barcodes-X3.tsv
sed -i "s/\([A-Z]\+\)\-1/\1\-X4/g" barcodes-X4.tsv
cat barcodes-X2.tsv barcodes-X3.tsv barcodes-X4.tsv > barcodes_merged.txt


