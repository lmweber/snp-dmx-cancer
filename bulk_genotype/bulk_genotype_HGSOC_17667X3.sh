#!/bin/bash

###############################################
# Shell script to genotype bulk RNA-seq samples
###############################################

# This script runs STAR and cellSNP (mode 2) to align and genotype matched bulk 
# RNA-seq samples from the same patients as our single-cell HGSOC samples. This 
# generates VCF files that can then be provided to cellSNP/Vireo to improve 
# demultiplexing performance.


# note: requires STAR index from previous script "create_STAR_index.sh"

# qsub -V -cwd -pe local 20 -l mem_free=10G,h_vmem=11G,h_fsize=100G bulk_genotype.sh


# --------
# Run STAR
# --------

# align reads

STAR \
--genomeDir ../../bulk/STAR_index \
--runThreadN 10 \
--readFilesIn ../../data/HGSOC/17667R/Fastq/17667X3_200214_A00421_0157_BHK52CDRXX_S41_L002_R1_001.fastq.gz ../../data/HGSOC/17667R/Fastq/17667X3_200214_A00421_0157_BHK52CDRXX_S41_L002_R2_001.fastq.gz \
--outFileNamePrefix ../../bulk/17667X3/STAR/ \
--readFilesCommand gunzip -c \
--outSAMtype BAM SortedByCoordinate \
--limitGenomeGenerateRAM 200000000000


# ---------
# Index BAM
# ---------

samtools index ../../bulk/17667X3/STAR/Aligned.sortedByCoord.out.bam


# --------------------
# Run cellSNP (mode 2)
# --------------------

# note: more stable to run interactively with "qrsh" instead of "qsub" on cluster

cellSNP \
-s ../../bulk/17667X3/STAR/Aligned.sortedByCoord.out.bam \
-O ../../bulk/17667X3/cellSNP \
-p 20 \
--minMAF=0.05 \
--minCOUNT=100 \
--UMItag=None

