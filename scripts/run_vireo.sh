#!/bin/bash

# -------------------------
# Shell script to run Vireo
# -------------------------

# run Vireo to demultiplex samples

# notes:
# - running Vireo using option 1 (without sample genotyping)
# - assuming known number of samples
# - requires cellSNP output from previous step

# for more details:
# - https://vireosnp.readthedocs.io/en/latest/manual.html

# runtime: ~10 min

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G run_vireo.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: output directory


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


vireo -c $4/cellSNP -N 3 -o $4/vireo --randSeed=123


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/vireo
echo runtime: $runtime seconds > $1/vireo/runtime_vireo.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/vireo
date > $2/vireo/timestamp_vireo.txt
# -----------------------------------

