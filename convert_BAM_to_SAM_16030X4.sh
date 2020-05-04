# Convert Cell Ranger output from BAM to SAM

# required to parse SAM output file


#qsub -V -cwd -l mem_free=100G,h_vmem=200G,h_fsize=300G convert_BAM_to_SAM.sh


samtools view -h 16030X4_HJTWLDMXX/outs/possorted_genome_bam.bam > 16030X4_HJTWLDMXX/outs/possorted_genome_bam.sam


