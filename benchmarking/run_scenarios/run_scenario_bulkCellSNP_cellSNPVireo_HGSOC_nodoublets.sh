#!/bin/bash

#######################################################
# Shell script to run simulation scenario (no doublets)
#######################################################

# These scripts run the selected demultiplexing tool (e.g. cellSNP/Vireo) for a 
# given simulation scenario (VCF file, dataset, percent doublets).

# Requires merged BAM file from benchmarking pipeline ("benchmarking/Snakefile").


# qsub -V -cwd -pe local 10 -l mem_free=5G,h_vmem=6G,h_fsize=100G run_scenario.sh


# start runtime
start=`date +%s`


# ------------------------------------------------------
# Scenario: bulk samples VCF from cellSNP, cellSNP/Vireo
# ------------------------------------------------------

# run cellSNP

mkdir -p ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkCellSNP_cellSNPVireo

# using recommended parameters for cellSNP
cellsnp-lite \
-s ../../../benchmarking/outputs/HGSOC/bam_merged/bam_merged.bam \
-b ../../../benchmarking/outputs/HGSOC/barcodes_merged/barcodes_merged.tsv \
-O ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkCellSNP_cellSNPVireo/cellSNP \
-R ../../../genotype/cellSNP_bulk_merged/cellSNP.cells-merged.vcf \
-p 10 \
--minMAF=0.1 \
--minCOUNT=20 \
--gzip


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkCellSNP_cellSNPVireo/runtimes
echo runtime: $runtime seconds > ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkCellSNP_cellSNPVireo/runtimes/runtime_bulkCellSNP_cellSNPVireo_cellSNP_HGSOC_nodoublets.txt

# start runtime
start=`date +%s`


# run Vireo

# note parameter for known number of samples (3 for HGSOC dataset, 6 for lung dataset)
vireo \
-c ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkCellSNP_cellSNPVireo/cellSNP \
-N 3 \
-o ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkCellSNP_cellSNPVireo/vireo \
--randSeed=123


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkCellSNP_cellSNPVireo/runtimes
echo runtime: $runtime seconds > ../../../benchmarking/scenarios/HGSOC/nodoublets/bulkCellSNP_cellSNPVireo/runtimes/runtime_bulkCellSNP_cellSNPVireo_vireo_HGSOC_nodoublets.txt

