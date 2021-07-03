#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l mem_free=5G,h_vmem=10G,h_fsize=100G


# ---------------------------
# Shell script to run cellSNP
# ---------------------------

# start runtime
start=`date +%s`

mkdir -p ../../supplementary_healthy/scenarios/nodoublets/1000GenomesFilt_cellSNPVireo

cellsnp-lite \
-s ../../supplementary_healthy/outputs/bam_merged/bam_merged.bam \
-b ../../supplementary_healthy/outputs/barcodes_merged/barcodes_merged.tsv \
-O ../../supplementary_healthy/scenarios/nodoublets/1000GenomesFilt_cellSNPVireo/cellSNP \
-R ../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf \
-p 10 \
--minMAF=0.1 \
--minCOUNT=20 \
--gzip


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../supplementary_healthy/runtimes/scenarios/nodoublets/cellSNP
echo runtime: $runtime seconds > ../../supplementary_healthy/runtimes/scenarios/nodoublets/cellSNP/runtime_cellSNP_1000GenomesFilt_cellSNPVireo.txt

