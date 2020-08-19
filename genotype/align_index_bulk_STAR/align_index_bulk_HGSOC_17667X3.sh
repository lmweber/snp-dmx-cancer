#!/bin/bash

#########################################################
# Shell script to align reads for bulk samples using STAR
#########################################################

# This script runs (i) STAR to align reads and (ii) samtools to index the 
# resulting BAM file for the bulk RNA-seq samples in our HGSOC dataset.


# note: requires STAR index from previous script "create_STAR_index.sh"

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=11G,h_fsize=100G align_bulk_STAR.sh


# start runtime
start=`date +%s`


# --------
# Run STAR
# --------

# align reads

STAR \
--genomeDir ../../../genotype/STAR_index \
--runThreadN 10 \
--readFilesIn ../../../data/HGSOC/17667R/Fastq/17667X3_200214_A00421_0157_BHK52CDRXX_S41_L002_R1_001.fastq.gz ../../../data/HGSOC/17667R/Fastq/17667X3_200214_A00421_0157_BHK52CDRXX_S41_L002_R2_001.fastq.gz \
--outFileNamePrefix ../../../genotype/17667X3/STAR/ \
--readFilesCommand gunzip -c \
--outSAMtype BAM SortedByCoordinate \
--limitGenomeGenerateRAM 200000000000


# ---------
# Index BAM
# ---------

samtools index ../../../genotype/17667X3/STAR/Aligned.sortedByCoord.out.bam


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../genotype/runtimes/align_index_bulk_STAR
echo runtime: $runtime seconds > ../../../genotype/runtimes/align_index_bulk_STAR/runtime_align_index_bulk_17667X3.txt

