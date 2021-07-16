#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=10G,h_fsize=100G


# -------------------------
# Shell script to run Vireo
# -------------------------

# start runtime
start=`date +%s`


# note parameter for known number of samples

vireo \
-c ../../../supplementary_healthy/scenarios/20pc/1000GenomesFilt_cellSNPVireo/cellSNP \
-N 5 \
-o ../../../supplementary_healthy/scenarios/20pc/1000GenomesFilt_cellSNPVireo/vireo \
--randSeed=123


# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
mkdir -p ../../../supplementary_healthy/runtimes/scenarios/20pc/vireo
echo runtime: $runtime seconds > ../../../supplementary_healthy/runtimes/scenarios/20pc/vireo/runtime_vireo_1000GenomesFilt_cellSNPVireo.txt

