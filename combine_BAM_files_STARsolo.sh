# Bash script to combine BAM files from STARsolo


# -----------------
# combine BAM files
# -----------------

#qsub -V -cwd -l mem_free=100G,h_vmem=110G,h_fsize=200G combine_BAM_files_STARsolo.sh


mkdir BAM_merged

samtools merge BAM_merged/BAM_merged_STARsolo.bam 16030X2/Aligned.out.bam 16030X3/Aligned.out.bam 16030X4/Aligned.out.bam


