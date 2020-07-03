# Bash script to run salmon alevin to create BAM files containing cell barcodes


# runtime: 1 hour

#qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=200G run_salmon_alevin.sh


# load required modules on JHPCE server
module load gcc/5.5.0
module load cmake/3.15.4

salmon alevin \
-l ISR \
-1 ../data/16030R/Fastq/16030X2_HJVMLDMXX/16030X2_HJVMLDMXX_S1_L002_R1_001.fastq.gz \
-2 ../data/16030R/Fastq/16030X2_HJVMLDMXX/16030X2_HJVMLDMXX_S1_L002_R2_001.fastq.gz \
--chromium \
-i ../data/index \
-p 10 \
-o ../outputs/16030X2_HJVMLDMXX/alevin_output \
--tgMap ../data/index/tx2gene.tsv \
--writeMappings=../outputs/16030X2_HJVMLDMXX/alevin_mappings/16030X2_HJVMLDMXX.sam \
--dumpMtx


