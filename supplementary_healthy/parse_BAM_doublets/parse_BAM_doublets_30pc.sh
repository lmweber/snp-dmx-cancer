#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=10G,h_fsize=200G


# ----------------------------------------------------------------
# Shell script to parse through merged BAM file to create doublets
# ----------------------------------------------------------------

# start runtime
start=`date +%s`


# parse through BAM file

# note hyphen for argument order
samtools view -h ../../../supplementary_healthy/outputs/bam_merged/bam_merged.bam | \
awk \
'NR==1 { next } FNR==NR { a[$1]=$2; next } (i=gensub(/.*CB\:Z\:([A-Za-z]+\-[A-Za-z0-9]+).*/, "\\1", 1, $0)) in a { gsub(i, a[i]) }1' \
../../../supplementary_healthy/scenarios/30pc/lookup_table_doublets_30pc.tsv - | \
samtools view -bo ../../../supplementary_healthy/scenarios/30pc/bam_merged_doublets_30pc.bam


# index BAM file

samtools index ../../../supplementary_healthy/scenarios/30pc/bam_merged_doublets_30pc.bam


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../supplementary_healthy/runtimes/parse_BAM_doublets
echo runtime: $runtime seconds > ../../../supplementary_healthy/runtimes/parse_BAM_doublets/runtime_parse_BAM_HGSOC_30pc.txt

