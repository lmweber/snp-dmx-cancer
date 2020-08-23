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


# ----------------------------------------------------------
# Scenario: 1000 Genomes Project VCF filtered, cellSNP/Vireo
# ----------------------------------------------------------

# run cellSNP

mkdir -p ../../../scenarios/doublets/HGSOC/20pc/1000GenomesFilt_cellSNPVireo

# using recommended parameters for cellSNP
cellSNP \
-s ../../../scenarios/doublets/HGSOC/20pc/bam_merged_doublets_HGSOC_20pc.bam \
-b ../../../scenarios/doublets/HGSOC/20pc/barcodes_merged_HGSOC_20pc.tsv \
-O ../../../scenarios/doublets/HGSOC/20pc/1000GenomesFilt_cellSNPVireo/cellSNP \
-R ../../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf \
-p 10 \
--minMAF=0.1 \
--minCOUNT=20


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../scenarios/doublets/HGSOC/20pc/1000GenomesFilt_cellSNPVireo/runtimes
echo runtime: $runtime seconds > ../../../scenarios/doublets/HGSOC/20pc/1000GenomesFilt_cellSNPVireo/runtimes/runtime_1000GenomesFilt_cellSNPVireo_cellSNP_HGSOC_20pc.txt

# start runtime
start=`date +%s`


# run Vireo

# note parameter for known number of samples (3 for HGSOC dataset, 6 for lung dataset)
vireo \
-c ../../../scenarios/doublets/HGSOC/20pc/1000GenomesFilt_cellSNPVireo/cellSNP \
-N 3 \
-o ../../../scenarios/doublets/HGSOC/20pc/1000GenomesFilt_cellSNPVireo/vireo \
--randSeed=123


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../scenarios/doublets/HGSOC/20pc/1000GenomesFilt_cellSNPVireo/runtimes
echo runtime: $runtime seconds > ../../../scenarios/doublets/HGSOC/20pc/1000GenomesFilt_cellSNPVireo/runtimes/runtime_1000GenomesFilt_cellSNPVireo_vireo_HGSOC_20pc.txt

