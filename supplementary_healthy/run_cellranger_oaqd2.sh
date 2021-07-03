#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l mem_free=5G,h_vmem=10G,h_fsize=200G


# -------------------------------
# Shell script to run Cell Ranger
# -------------------------------

# start runtime
start=`date +%s`

# set working directory for output path
cwd=$(pwd)
mkdir -p ../../supplementary_healthy/outputs/oaqd2
cd ../../supplementary_healthy/outputs/oaqd2


cellranger count --id=oaqd2 \
--description=oaqd2 \
--transcriptome=../../../data/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=../../../data/souporcell/fastq/oaqd2 \
--sample=oaqd2 \
--nosecondary \
--jobmode=local \
--localcores=10 \
--localmem=50 \
--localvmem=100


# restore working directory
cd $cwd

# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../supplementary_healthy/runtimes/cellranger
echo runtime: $runtime seconds > ../../supplementary_healthy/runtimes/cellranger/runtime_cellranger_oaqd2.txt

