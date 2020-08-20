#!/bin/bash

######################################################################
# Shell script to genotype bulk RNA-seq samples using bcftools mpileup
######################################################################

# This script runs bcftools mpileup to genotype and generate VCF files for the 
# bulk RNA-seq samples in our HGSOC dataset.


# note: requires BAM files from previous scripts "align_index_bulk_HGSOC_17667XX.sh"

# qsub -V -cwd -l mem_free=2G,h_vmem=3G,h_fsize=100G genotype_bulk_HGSOC_bcftools.sh


# start runtime
start=`date +%s`


# --------------------
# Run bcftools mpileup
# --------------------

mkdir -p ../../../genotype/bcftools

bcftools mpileup -Ou \
-f ../../../data/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
../../../genotype/17667X1/STAR/Aligned.sortedByCoord.out.bam \
../../../genotype/17667X2/STAR/Aligned.sortedByCoord.out.bam \
../../../genotype/17667X3/STAR/Aligned.sortedByCoord.out.bam | \
bcftools call -mv -Ov \
-o ../../../genotype/bcftools/bcftools_HGSOC.vcf


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../genotype/runtimes/genotype_bulk_bcftools
echo runtime: $runtime seconds > ../../../genotype/runtimes/genotype_bulk_bcftools/runtime_genotype_bulk_HGSOC_bcftools.txt

