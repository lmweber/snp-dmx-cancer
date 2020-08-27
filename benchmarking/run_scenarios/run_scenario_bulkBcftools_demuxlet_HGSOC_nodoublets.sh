#!/bin/bash

#######################################################
# Shell script to run simulation scenario (no doublets)
#######################################################

# These scripts run the selected demultiplexing tool (e.g. cellSNP/Vireo) for a 
# given simulation scenario (VCF file, dataset, percent doublets).

# Requires merged BAM file from benchmarking pipeline ("benchmarking/Snakefile").


# qsub -V -cwd -l mem_free=5G,h_vmem=6G,h_fsize=100G run_scenario.sh


# start runtime
start=`date +%s`


# --------------------------------------------------
# Scenario: bulk samples VCF from bcftools, demuxlet
# --------------------------------------------------

mkdir -p ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkBcftools_demuxlet

# run demuxlet

demuxlet \
--sam ../../../benchmarking/scenarios/outputs/HGSOC/bam_merged.bam \
--group-list ../../../benchmarking/scenarios/outputs/HGSOC/barcodes_merged.tsv \
--alpha 0 --alpha 0.5 \
--vcf ../../../genotype/bcftools/bcftools_HGSOC_rehead.vcf \
--field GT \
--out ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkBcftools_demuxlet/demuxlet


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkBcftools_demuxlet/runtimes
echo runtime: $runtime seconds > ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkBcftools_demuxlet/runtimes/runtime_bulkBcftools_demuxlet_HGSOC_nodoublets.txt

