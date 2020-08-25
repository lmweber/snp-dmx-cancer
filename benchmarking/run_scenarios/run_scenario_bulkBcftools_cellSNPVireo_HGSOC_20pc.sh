#!/bin/bash

##################################################
# Shell script to run doublets simulation scenario
##################################################

# These scripts run the selected demultiplexing tool (e.g. cellSNP/Vireo) for a 
# given doublets simulation scenario (VCF file, dataset, percent doublets).

# Requires the modified BAM file from the previous script "parse_BAM_doublets_X.sh".


# qsub -V -cwd -pe local 10 -l mem_free=5G,h_vmem=6G,h_fsize=100G run_scenario.sh


# start runtime
start=`date +%s`


# -------------------------------------------------------
# Scenario: bulk samples VCF from bcftools, cellSNP/Vireo
# -------------------------------------------------------

# run cellSNP

mkdir -p ../../../benchmarking/scenarios/HGSOC/20pc/bulkBcftools_cellSNPVireo

# using recommended parameters for cellSNP
cellSNP \
-s ../../../benchmarking/scenarios/HGSOC/20pc/bam_merged_doublets_HGSOC_20pc.bam \
-b ../../../benchmarking/scenarios/HGSOC/20pc/barcodes_merged_HGSOC_20pc.tsv \
-O ../../../benchmarking/scenarios/HGSOC/20pc/bulkBcftools_cellSNPVireo/cellSNP \
-R ../../../genotype/bcftools/bcftools_HGSOC_rehead.vcf \
-p 10 \
--minMAF=0.1 \
--minCOUNT=20


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../benchmarking/scenarios/HGSOC/20pc/bulkBcftools_cellSNPVireo/runtimes
echo runtime: $runtime seconds > ../../../benchmarking/scenarios/HGSOC/20pc/bulkBcftools_cellSNPVireo/runtimes/runtime_bulkBcftools_cellSNPVireo_cellSNP_HGSOC_20pc.txt

# start runtime
start=`date +%s`


# run Vireo

# note parameter for known number of samples (3 for HGSOC dataset, 6 for lung dataset)
vireo \
-c ../../../benchmarking/scenarios/HGSOC/20pc/bulkBcftools_cellSNPVireo/cellSNP \
-N 3 \
-o ../../../benchmarking/scenarios/HGSOC/20pc/bulkBcftools_cellSNPVireo/vireo \
--randSeed=123


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../benchmarking/scenarios/HGSOC/20pc/bulkBcftools_cellSNPVireo/runtimes
echo runtime: $runtime seconds > ../../../benchmarking/scenarios/HGSOC/20pc/bulkBcftools_cellSNPVireo/runtimes/runtime_bulkBcftools_cellSNPVireo_vireo_HGSOC_20pc.txt

