# Bash script to run STARsolo to create BAM files for scRNA-seq samples

# STARsolo is a faster and less buggy alternative to Cell Ranger
# see https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

# requires STAR genome index from previous script


# ------------
# run STARsolo
# ------------

#qsub -V -cwd -pe local 20 -l mem_free=10G,h_vmem=11G,h_fsize=100G run_STARsolo.sh

# note: R2 is before R1 in --readFilesIn


# STAR --genomeDir ../data/STAR_index \
# --readFilesCommand zcat \
# --outFileNamePrefix ./16030X2/ \
# --runThreadN 20 \
# --readFilesIn ../data/16030R/Fastq/16030X2_HJVMLDMXX/16030X2_HJVMLDMXX_S1_L002_R2_001.fastq.gz ../data/16030R/Fastq/16030X2_HJVMLDMXX/16030X2_HJVMLDMXX_S1_L002_R1_001.fastq.gz \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ../data/cellranger_whitelist/3M-february-2018.txt \
# --soloUMIlen 12 \
# --soloUMIfiltering MultiGeneUMI \
# --soloCBmatchWLtype 1MM_multi_pseudocounts

# STAR --genomeDir ../data/STAR_index \
# --readFilesCommand zcat \
# --outFileNamePrefix ./16030X3/ \
# --runThreadN 20 \
# --readFilesIn ../data/16030R/Fastq/16030X3_HJTWLDMXX/16030X3_HJTWLDMXX_S1_L001_R2_001.fastq.gz ../data/16030R/Fastq/16030X3_HJTWLDMXX/16030X3_HJTWLDMXX_S1_L001_R1_001.fastq.gz  \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ../data/cellranger_whitelist/3M-february-2018.txt \
# --soloUMIlen 12 \
# --soloUMIfiltering MultiGeneUMI \
# --soloCBmatchWLtype 1MM_multi_pseudocounts

STAR --genomeDir ../data/STAR_index \
--readFilesCommand zcat \
--outFileNamePrefix ./16030X4/ \
--runThreadN 20 \
--readFilesIn ../data/16030R/Fastq/16030X4_HJTWLDMXX/16030X4_HJTWLDMXX_S2_L001_R2_001.fastq.gz ../data/16030R/Fastq/16030X4_HJTWLDMXX/16030X4_HJTWLDMXX_S2_L001_R1_001.fastq.gz  \
--soloType CB_UMI_Simple \
--soloCBwhitelist ../data/cellranger_whitelist/3M-february-2018.txt \
--soloUMIlen 12 \
--soloUMIfiltering MultiGeneUMI \
--soloCBmatchWLtype 1MM_multi_pseudocounts


