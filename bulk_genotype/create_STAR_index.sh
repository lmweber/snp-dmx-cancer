#!/bin/bash

###################################
# Shell script to create STAR index
###################################

# This script runs STAR to create genome index files, required by the following 
# scripts "bulk_genotype_HGSOC_17667XX.sh".


# note: creating STAR index requires up to 200GB of memory

# qsub -V -cwd -pe local 20 -l mem_free=10G,h_vmem=11G,h_fsize=200G create_STAR_index.sh


# -----------------
# Create STAR index
# -----------------

# create genome index

STAR \
--runMode genomeGenerate \
--runThreadN 20 \
--genomeDir ../../bulk/STAR_index \
--genomeFastaFiles ../../data/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
--sjdbGTFfile ../../data/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
--limitGenomeGenerateRAM 200000000000

