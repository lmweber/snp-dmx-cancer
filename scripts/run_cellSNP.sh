#!/bin/bash

# ---------------------------
# Shell script to run cellSNP
# ---------------------------

# run cellSNP to genotype cells

# notes:
# - requires merged BAM file and merged cell barcodes file, each containing unique sample IDs
# - running cellSNP in mode 1 (using .vcf file containing common variants)

# for more details:
# - https://vireosnp.readthedocs.io/en/latest/genotype.html
# - https://github.com/single-cell-genetics/cellSNP

# runtime: ~30 min

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G run_cellSNP.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: data directory
# $5: output directory


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------

# note: vcf file needs to be uncompressed
# if still in "vcf.bgz" format then uncompress as follows:
# mv file.vcf.bgz file.vcf.gz
# gunzip file.vcf.gz


cellSNP \
-s $5/bam_merged/bam_merged.bam \
-b $5/barcodes_merged/barcodes_merged.txt \
-O $5/cellSNP \
-R $4/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs_transcriptomic_v2.vcf \
-p $3 \
--minMAF=0.05


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/cellSNP
echo runtime: $runtime seconds > $1/cellSNP/runtime_cellSNP.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/cellSNP
date > $2/cellSNP/timestamp_cellSNP.txt
# -----------------------------------

