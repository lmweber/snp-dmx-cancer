# --------------------------------
# Bash script to run salmon alevin
# --------------------------------

# runtime: ~1 hour

# qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=200G run_salmon_alevin.sh

# arguments:
# $1: sample ID
# $2: directory for runtimes
# $3: directory for timestamp files
# $4: number of threads
# $5: filename for FASTQ file 1
# $6: filename for FASTQ file 2
# $7: salmon index directory
# $8: output directory


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


# load required modules on our JHPCE cluster
module load gcc/5.5.0
module load cmake/3.15.4

# note: run salmon alevin with option '--writeMappings' to output SAM file containing cell barcodes

salmon alevin \
-l ISR \
-1 $5 \
-2 $6 \
--chromium \
-i $7 \
-p $4 \
-o $8/alevin_output \
--tgMap $7/tx2gene.tsv \
--writeMappings=$8/alevin_mappings/$1.sam


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p $2/HGSOC/alevin
echo runtime: $runtime seconds > $2/HGSOC/alevin/runtime_alevin_$1.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
mkdir -p $3/HGSOC/alevin
date > $3/HGSOC/alevin/timestamp_alevin_$1.txt
# -----------------------------------

