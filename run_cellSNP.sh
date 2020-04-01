# Bash script to run cellSNP to genotype cells

# see https://vireosnp.readthedocs.io/en/latest/genotype.html
# and https://github.com/single-cell-genetics/cellSNP


# ------------
# installation
# ------------

# using Python module on JHPCE

#module load python/3.7.3

#pip install cellSNP --user

#pip list -v


# -------------
# cell barcodes
# -------------

# get cell barcodes from Cell Ranger outputs from each sample, and combine into a single file

# cd CellRanger/barcodes
# cp ../16030X2_HJVMLDMXX/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ./barcodes-16030X2_HJVMLDMXX.tsv.gz
# cp ../16030X3_HJTWLDMXX/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ./barcodes-16030X3_HJTWLDMXX.tsv.gz
# cp ../16030X4_HJTWLDMXX/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ./barcodes-16030X4_HJTWLDMXX.tsv.gz
# gunzip *
# cat barcodes-16030X2_HJVMLDMXX.tsv barcodes-16030X3_HJTWLDMXX.tsv barcodes-16030X4_HJTWLDMXX.tsv > barcodes_merged.txt
# # remove '-1' from end of barcodes
# sed 's/-[0-9]*$//' barcodes_merged.txt > barcodes_merged_clean.txt


# -----------
# run cellSNP
# -----------

# "mode 1": requires a single merged BAM file from multiple scRNA-seq samples


#cd /fastscratch/myscratch/lweber/cellSNP

#qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=200G run_cellSNP.sh


BAM=../CellRanger/BAM_merged/BAM_merged_CellRanger.bam
#BARCODE=../CellRanger/barcodes/barcodes_merged_clean.txt
BARCODE=../data/cellranger_whitelist/3M-february-2018.txt
REGION_VCF=../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf


module load python/3.7.3

cellSNP -s $BAM -b $BARCODE -O $OUT_DIR -R $REGION_VCF -p 20 --minMAF 0.1 --minCOUNT 20


