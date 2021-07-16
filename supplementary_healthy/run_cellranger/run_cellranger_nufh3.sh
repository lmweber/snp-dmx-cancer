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
mkdir -p ../../../supplementary_healthy/outputs/nufh3
cd ../../../supplementary_healthy/outputs/nufh3


cellranger count --id=nufh3 \
--description=nufh3 \
--transcriptome=../../../../data/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=../../../../data/souporcell/fastq/nufh3 \
--sample=nufh3 \
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
mkdir -p ../../../supplementary_healthy/runtimes/cellranger
echo runtime: $runtime seconds > ../../../supplementary_healthy/runtimes/cellranger/runtime_cellranger_nufh3.txt

