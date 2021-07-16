#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=10G,h_fsize=500G


# ------------------------------------------------
# Shell script to merge and index parsed BAM files
# ------------------------------------------------

# start runtime
start=`date +%s`

mkdir -p ../../../supplementary_healthy/outputs/bam_merged


samtools merge ../../../supplementary_healthy/outputs/bam_merged/bam_merged.bam \
../../../supplementary_healthy/outputs/euts1/euts1/outs/possorted_genome_bam_parsed.bam \
../../../supplementary_healthy/outputs/nufh3/nufh3/outs/possorted_genome_bam_parsed.bam \
../../../supplementary_healthy/outputs/babz3/babz3/outs/possorted_genome_bam_parsed.bam \
../../../supplementary_healthy/outputs/oaqd2/oaqd2/outs/possorted_genome_bam_parsed.bam \
../../../supplementary_healthy/outputs/ieki3/ieki3/outs/possorted_genome_bam_parsed.bam


samtools index ../../../supplementary_healthy/outputs/bam_merged/bam_merged.bam


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../supplementary_healthy/runtimes/merge_and_index_BAM
echo runtime: $runtime seconds > ../../../supplementary_healthy/runtimes/merge_and_index_BAM/runtime_merge_and_index_BAM.txt

