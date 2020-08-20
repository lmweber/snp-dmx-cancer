#!/bin/bash

################################################################
# Shell script to update samples names in VCF file from bcftools
################################################################

# This script runs bcftools reheader to update sample names in the VCF file 
# previously generated with bcftools mpileup.


# note: requires VCF file from previous script "genotype_bulk_HGSOC_bcftools.sh"

# qsub -V -cwd -l mem_free=2G,h_vmem=3G,h_fsize=100G genotype_bulk_HGSOC_bcftools_reheader.sh


# start runtime
start=`date +%s`


# ---------------------
# Run bcftools reheader
# ---------------------

# create file containing updated sample names
echo "../../../genotype/17667X1/STAR/Aligned.sortedByCoord.out.bam 17667X1" > ../../../genotype/bcftools/sample_names_bulk.txt
echo "../../../genotype/17667X2/STAR/Aligned.sortedByCoord.out.bam 17667X2" >> ../../../genotype/bcftools/sample_names_bulk.txt
echo "../../../genotype/17667X3/STAR/Aligned.sortedByCoord.out.bam 17667X3" >> ../../../genotype/bcftools/sample_names_bulk.txt

# update sample names in VCF file
cp ../../../genotype/bcftools/bcftools_HGSOC.vcf ../../../genotype/bcftools/bcftools_HGSOC_bak.vcf
bcftools reheader -s sample_names_bulk.txt bcftools_HGSOC_bak.vcf > bcftools_HGSOC.vcf


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../genotype/runtimes/genotype_bulk_bcftools
echo runtime: $runtime seconds > ../../../genotype/runtimes/genotype_bulk_bcftools/runtime_genotype_bulk_HGSOC_bcftools_reheader.txt

