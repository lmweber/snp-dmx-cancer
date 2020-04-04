# Bash script to merge BAM files from Cell Ranger using 'synth_pool.py' script from Vireo authors


# -------------------------------
# unzip Cell Ranger barcode files
# -------------------------------

# cd 16030X2_HJVMLDMXX/outs/filtered_feature_bc_matrix
# gunzip -c barcodes.tsv.gz > barcodes.tsv
# 
# cd 16030X3_HJTWLDMXX/outs/filtered_feature_bc_matrix
# gunzip -c barcodes.tsv.gz > barcodes.tsv
# 
# cd 16030X4_HJTWLDMXX/outs/filtered_feature_bc_matrix
# gunzip -c barcodes.tsv.gz > barcodes.tsv


# -------------------------------------
# merge BAM files using 'synth_pool.py'
# -------------------------------------

#qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=200G merge_BAM_files_synth_pool.sh

mkdir BAM_merged

python ../synth_pool/synth_pool.py \
-s 16030X2_HJVMLDMXX/outs/possorted_genome_bam.bam,16030X3_HJTWLDMXX/outs/possorted_genome_bam.bam,16030X4_HJTWLDMXX/outs/possorted_genome_bam.bam \
-b 16030X2_HJVMLDMXX/outs/filtered_feature_bc_matrix/barcodes.tsv,16030X3_HJTWLDMXX/outs/filtered_feature_bc_matrix/barcodes.tsv,16030X4_HJTWLDMXX/outs/filtered_feature_bc_matrix/barcodes.tsv \
-o BAM_merged \
-p 10


