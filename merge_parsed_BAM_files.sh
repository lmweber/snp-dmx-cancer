# Merge parsed BAM files from Cell Ranger


# notes:
# assumes BAM files have been parsed to contain unique sample IDs in cell barcodes


#qsub -V -cwd -l mem_free=100G,h_vmem=200G,h_fsize=200G merge_BAM_files.sh


mkdir bam_merged

samtools merge bam_merged/bam_merged.bam 16030X2_HJVMLDMXX/outs/possorted_genome_bam_parsed.bam 16030X3_HJTWLDMXX/outs/possorted_genome_bam_parsed.bam 16030X4_HJTWLDMXX/outs/possorted_genome_bam_parsed.bam


