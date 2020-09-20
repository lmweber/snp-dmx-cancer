#!/bin/bash

# -------------------------
# Shell script to run Vireo
# -------------------------

# run Vireo to demultiplex samples

# notes:
# - assuming known number of samples
# - requires cellSNP output from previous step

# for more details:
# - https://vireosnp.readthedocs.io/en/latest/manual.html

# runtime: ~30 min

# qsub -V -cwd -l mem_free=5G,h_vmem=10G,h_fsize=100G run_vireo.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: output directory
# $5: number of samples
# $6: dataset name for simulation scenario
# $7: percentage of doublets for simulation scenario (formatted as e.g. "20pc")


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------

# note parameter for known number of samples (3 for HGSOC dataset)

vireo -c $4/$6/doublets_sims/$7/cellSNP -N $5 -o $4/$6/doublets_sims/$7/vireo --randSeed=123


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $1/$6/doublets_sims/$7/vireo
echo runtime: $runtime seconds > $1/$6/doublets_sims/$7/vireo/runtime_vireo_$6_$7.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $2/$6/doublets_sims/$7/vireo
date > $2/$6/doublets_sims/$7/vireo/timestamp_vireo_$6_$7.txt
# -----------------------------------

