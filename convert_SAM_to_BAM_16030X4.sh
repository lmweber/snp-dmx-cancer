# Bash script to convert STARsolo output from SAM to BAM files

# requires STARsolo outputs from previous scripts


# ------------------
# convert SAM to BAM
# ------------------

#qsub -V -l mem_free=20G,h_vmem=21G,h_fsize=200G convert_SAM_to_BAM.sh


# samtools view -bS 16030X2/Aligned.out.sam > 16030X2/Aligned.out.bam

# samtools view -bS 16030X3/Aligned.out.sam > 16030X3/Aligned.out.bam

samtools view -bS 16030X4/Aligned.out.sam > 16030X4/Aligned.out.bam


