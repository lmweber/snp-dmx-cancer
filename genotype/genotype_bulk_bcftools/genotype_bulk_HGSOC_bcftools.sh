#!/bin/bash

######################################################################
# Shell script to genotype bulk RNA-seq samples using bcftools mpileup
######################################################################

# This script runs bcftools mpileup to genotype and generate VCF files for the 
# bulk RNA-seq samples in our HGSOC dataset.


# note: requires BAM files from previous scripts "align_bulk_HGSOC_17667XX.sh"

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=11G,h_fsize=100G genotype_bulk_HGSOC_bcftools.sh


# --------------------
# Run bcftools mpileup
# --------------------

bcftools mpileup -Ou \
-f ../../../data/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
../../../genotype/17667X1/STAR/Aligned.sortedByCoord.out.bam \
../../../genotype/17667X2/STAR/Aligned.sortedByCoord.out.bam \
../../../genotype/17667X3/STAR/Aligned.sortedByCoord.out.bam | \
bcftools call -mv -Ov \
-o ../../../genotype/bcftools/bcftools_HGSOC.vcf

