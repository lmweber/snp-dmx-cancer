# Bash script to combine BAM files from Cell Ranger


# -----------------
# combine BAM files
# -----------------

#qsub -V -cwd -l mem_free=100G,h_vmem=110G,h_fsize=200G combine_BAM_files_CellRanger.sh


mkdir BAM_merged

samtools merge BAM_merged/BAM_merged_CellRanger.bam 16030X2_HJVMLDMXX/outs/possorted_genome_bam.bam 16030X3_HJTWLDMXX/outs/possorted_genome_bam.bam 16030X4_HJTWLDMXX/outs/possorted_genome_bam.bam


