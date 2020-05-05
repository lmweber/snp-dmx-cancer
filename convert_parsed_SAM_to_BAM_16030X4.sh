# Convert parsed SAM file back to BAM format

# required otherwise merged SAM file is too large


#qsub -V -cwd -l mem_free=100G,h_vmem=300G,h_fsize=300G convert_parsed_SAM_to_BAM.sh


samtools view -S -b 16030X4_HJTWLDMXX/outs/possorted_genome_bam.sam > 16030X4_HJTWLDMXX/outs/possorted_genome_bam_parsed.bam


