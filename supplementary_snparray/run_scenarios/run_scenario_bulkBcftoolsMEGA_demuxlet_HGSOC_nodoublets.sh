#!/bin/bash

##################################################
# Shell script to run doublets simulation scenario
##################################################

# These scripts run the selected demultiplexing tool (e.g. cellSNP/Vireo) for a 
# given doublets simulation scenario (VCF file, dataset, percent doublets).

# Requires the modified BAM file from the previous script "parse_BAM_doublets_X.sh".


# qsub -l mem_free=5G,h_vmem=6G,h_fsize=100G run_scenario.sh


# start runtime
start=`date +%s`


# --------------------------------------------------------------------------
# Scenario: bulk samples VCF from bcftools (subset MEGA SNP array), demuxlet
# --------------------------------------------------------------------------

mkdir -p ../../../supplementary_snparray/scenarios/HGSOC/nodoublets/bulkBcftoolsMEGA_demuxlet

# run demuxlet

demuxlet \
--sam ../../../benchmarking/outputs/HGSOC/bam_merged/bam_merged.bam \
--group-list ../../../benchmarking/outputs/HGSOC/barcodes_merged/barcodes_merged.tsv \
--alpha 0 --alpha 0.5 \
--vcf ../../../genotype/MEGA/subset_MEGA_bcftools.vcf \
--field GT \
--out ../../../supplementary_snparray/scenarios/HGSOC/nodoublets/bulkBcftoolsMEGA_demuxlet/demuxlet


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../supplementary_snparray/scenarios/HGSOC/nodoublets/bulkBcftoolsMEGA_demuxlet/runtimes
echo runtime: $runtime seconds > ../../../supplementary_snparray/scenarios/HGSOC/nodoublets/bulkBcftoolsMEGA_demuxlet/runtimes/runtime_bulkBcftoolsMEGA_demuxlet_cellSNP_HGSOC_nodoublets.txt

