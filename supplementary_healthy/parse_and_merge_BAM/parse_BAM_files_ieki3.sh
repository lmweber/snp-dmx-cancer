#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=10G,h_fsize=200G


# ------------------------------------------------------------------
# Shell script to parse BAM files to add sample IDs to cell barcodes
# ------------------------------------------------------------------

# start runtime
start=`date +%s`


samtools view -h ../../../supplementary_healthy/outputs/ieki3/ieki3/outs/possorted_genome_bam.bam | \
sed "s|\(CB\:Z\:[A-Z]\+\)\-1|\1\-ieki3|g" | \
samtools view -bo ../../../supplementary_healthy/outputs/ieki3/ieki3/outs/possorted_genome_bam_parsed.bam


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../supplementary_healthy/runtimes/parse_BAM_files
echo runtime: $runtime seconds > ../../../supplementary_healthy/runtimes/parse_BAM_files/runtime_parse_BAM_files_ieki3.txt

