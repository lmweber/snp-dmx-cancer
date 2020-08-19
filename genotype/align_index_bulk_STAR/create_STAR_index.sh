#!/bin/bash

###################################
# Shell script to create STAR index
###################################

# This script runs STAR to create genome index files, required by the following 
# scripts to run STAR to align reads for the bulk RNA-seq samples.


# note: creating STAR index requires up to 200GB of memory

# qsub -V -cwd -pe local 10 -l mem_free=20G,h_vmem=21G,h_fsize=100G create_STAR_index.sh


# start runtime
start=`date +%s`


# -----------------
# Create STAR index
# -----------------

# create genome index

STAR \
--runMode genomeGenerate \
--runThreadN 10 \
--genomeDir ../../../genotype/STAR_index \
--genomeFastaFiles ../../../data/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
--sjdbGTFfile ../../../data/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
--limitGenomeGenerateRAM 200000000000


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../genotype/runtimes/align_index_bulk_STAR
echo runtime: $runtime seconds > ../../../genotype/runtimes/align_index_bulk_STAR/runtime_create_STAR_index.txt

