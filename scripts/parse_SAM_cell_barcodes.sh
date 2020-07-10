#!/bin/bash

# ------------------------------------------------------------------
# Shell script to parse SAM files to add sample IDs to cell barcodes
# ------------------------------------------------------------------

# notes:
# - syntax to search and replace: sed -i "s/regexp/replacement/g"
# - using double quotes allows variables inside sed expression
# - using alternative separator | allows slashes in variable inside sed expression
# - regular expression matches cell barcode with format "CB:Z:AGCTTCCAGCAGTCTT"
# - we add a suffix with the sample ID, e.g. "-X2" for sample X2

# runtime: ~2 hours

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G parse_SAM_cell_barcodes.sh

# arguments:
# $1: sample ID
# $2: directory for runtimes
# $3: directory for timestamp files
# $4: number of threads
# $5: short sample ID
# $6: input/output directory


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


sed -i "s|\(CB\:Z\:[A-Z]\+\)|\1\-$5|g" $6/$1.sam


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $2/HGSOC/parse_SAM_barcodes
echo runtime: $runtime seconds > $2/HGSOC/parse_SAM_barcodes/runtime_parse_SAM_barcodes_$1.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $3/HGSOC/parse_SAM_barcodes
date > $3/HGSOC/parse_SAM_barcodes/timestamp_parse_SAM_barcodes_$1.txt
# -----------------------------------

