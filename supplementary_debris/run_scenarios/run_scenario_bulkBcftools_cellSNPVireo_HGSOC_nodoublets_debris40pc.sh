#!/bin/bash

################################################
# Shell script to run debris simulation scenario
################################################

# These scripts run the selected demultiplexing tool (e.g. cellSNP/Vireo) for a 
# given debris simulation scenario (VCF file, dataset, percent doublets, percent 
# debris).

# Requires the modified BAM file from the previous script "parse_BAM_X_debris.sh".


# qsub -cwd -pe local 4 -l mem_free=5G,h_vmem=6G,h_fsize=100G run_scenario.sh


# start runtime
start=`date +%s`


# -------------------------------------------------------
# Scenario: bulk samples VCF from bcftools, cellSNP/Vireo
# -------------------------------------------------------

# run cellSNP

mkdir -p ../../../supplementary_debris/scenarios/HGSOC/nodoublets/debris40pc/bulkBcftools_cellSNPVireo

# using recommended parameters for cellSNP
cellsnp-lite \
-s ../../../supplementary_debris/scenarios/HGSOC/nodoublets/bam_merged_HGSOC_nodoublets_debris40pc.bam \
-b ../../../supplementary_debris/scenarios/HGSOC/nodoublets/barcodes_merged_HGSOC_nodoublets_debris40pc.tsv \
-O ../../../supplementary_debris/scenarios/HGSOC/nodoublets/debris40pc/bulkBcftools_cellSNPVireo/cellSNP \
-R ../../../genotype/bcftools/bcftools_HGSOC_rehead.vcf \
-p 4 \
--minMAF=0.1 \
--minCOUNT=20 \
--gzip


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../supplementary_debris/scenarios/HGSOC/nodoublets/debris40pc/bulkBcftools_cellSNPVireo/runtimes
echo runtime: $runtime seconds > ../../../supplementary_debris/scenarios/HGSOC/nodoublets/debris40pc/bulkBcftools_cellSNPVireo/runtimes/runtime_bulkBcftools_cellSNPVireo_cellSNP_HGSOC_nodoublets_debris40pc.txt

# start runtime
start=`date +%s`


# run Vireo

# note parameter for known number of samples (3 for HGSOC dataset, 6 for lung dataset)
vireo \
-c ../../../supplementary_debris/scenarios/HGSOC/nodoublets/debris40pc/bulkBcftools_cellSNPVireo/cellSNP \
-N 3 \
-o ../../../supplementary_debris/scenarios/HGSOC/nodoublets/debris40pc/bulkBcftools_cellSNPVireo/vireo \
--randSeed=123


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../supplementary_debris/scenarios/HGSOC/nodoublets/debris40pc/bulkBcftools_cellSNPVireo/runtimes
echo runtime: $runtime seconds > ../../../supplementary_debris/scenarios/HGSOC/nodoublets/debris40pc/bulkBcftools_cellSNPVireo/runtimes/runtime_bulkBcftools_cellSNPVireo_vireo_HGSOC_nodoublets_debris40pc.txt

