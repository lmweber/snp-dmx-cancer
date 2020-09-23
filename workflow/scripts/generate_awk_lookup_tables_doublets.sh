#!/bin/bash

# -----------------------------------------------------------------------------------
# Shell wrapper for R script to generate awk lookup tables and updated barcodes files
# -----------------------------------------------------------------------------------

# Shell wrapper script to run R script "generate_awk_lookup_tables_doublets.R" 
# on computing cluster, required on our cluster since R module cannot be loaded 
# directly in Snakemake rule. If R is installed locally, this shell script is 
# not required, and the R script can be run directly instead.

# runtime: seconds

# qsub -V -cwd -l mem_free=2G,h_vmem=3G,h_fsize=100G generate_awk_lookup_tables_doublets.sh

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: output directory


module load conda_R/4.0

Rscript generate_awk_lookup_tables_doublets.R $1 $2 $3

