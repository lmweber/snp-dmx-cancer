# Index merged BAM file

# index merged BAM file that was created in previous steps


# runtime: ~1 hour

#qsub -V -cwd -l mem_free=100G,h_vmem=300G,h_fsize=300G index_merged_BAM.sh


samtools index bam_merged/bam_merged.bam


