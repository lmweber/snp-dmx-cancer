#!/bin/bash

##################################################
# Shell script to run doublets simulation scenario
##################################################

# These scripts run the selected demultiplexing tool (e.g. cellSNP/Vireo) for a 
# given doublets simulation scenario (VCF file, dataset, percent doublets).

# Requires the modified BAM file from the previous script "parse_BAM_doublets_X.sh".


# qsub -V -cwd -l mem_free=5G,h_vmem=6G,h_fsize=100G run_scenario.sh


# start runtime
start=`date +%s`


# --------------------------------------------------
# Scenario: bulk samples VCF from bcftools, demuxlet
# --------------------------------------------------

mkdir -p ../../../benchmarking/scenarios/HGSOC/30pc/bulkBcftools_demuxlet

# run demuxlet

demuxlet \
--sam ../../../benchmarking/scenarios/HGSOC/30pc/bam_merged_doublets_HGSOC_30pc.bam \
--group-list ../../../benchmarking/scenarios/HGSOC/30pc/barcodes_merged_HGSOC_30pc.tsv \
--alpha 0 --alpha 0.5 \
--vcf ../../../genotype/bcftools/bcftools_HGSOC_rehead.vcf \
--field GT \
--out ../../../benchmarking/scenarios/HGSOC/30pc/bulkBcftools_demuxlet/demuxlet


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../benchmarking/scenarios/HGSOC/30pc/bulkBcftools_demuxlet/runtimes
echo runtime: $runtime seconds > ../../../benchmarking/scenarios/HGSOC/30pc/bulkBcftools_demuxlet/runtimes/runtime_bulkBcftools_demuxlet_HGSOC_30pc.txt

