#!/bin/bash

# --------------------------------
# Shell script to run salmon index
# --------------------------------

# runtime: ~5 min

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G run_salmon_index.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: sample ID
# $5: salmon index directory


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


# (1) download files (if required)

wget -P $5 ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.pc_transcripts.fa.gz
wget -P $5 ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz


# (2) run Python scripts (from Ariel Hippen Anderson)

# note: assumes valid python3 installation in path, with some required packages e.g. pandas

# make gene by transcript mapping file
python3 scripts/make_tx2gene.py -d $5 -f gencode.v34.pc_transcripts.fa.gz -g gencode.v34.annotation.gtf.gz

# run make_gene2symbol.py
python3 scripts/make_gene2symbol.py -d $5 -f gencode.v34.pc_transcripts.fa.gz -g gencode.v34.annotation.gtf.gz


# (3) run salmon index

# load required modules on our JHPCE cluster
module load gcc/5.5.0
module load cmake/3.15.4

salmon index -i $5 --gencode -p $3 -t $5/gencode.v34.pc_transcripts.fa.gz


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/salmon_index
echo runtime: $runtime seconds > $1/salmon_index/runtime_salmon_index.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/salmon_index
date > $2/salmon_index/timestamp_salmon_index.txt
# -----------------------------------

