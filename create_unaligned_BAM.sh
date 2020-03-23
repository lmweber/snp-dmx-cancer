# Bash script to create unaligned BAM files from FASTQ files
# (Vireo requires unaligned BAM format)

# following this tutorial:
# https://gatkforums.broadinstitute.org/gatk/discussion/6484/how-to-generate-an-unmapped-bam-from-fastq-or-aligned-bam


# --------
# tutorial
# --------

# cd ../FastqToSam
# 
# java -jar ../picard/picard.jar FastqToSam \
#     FASTQ=tutorial/tutorial_6484_FastqToSam/6484_snippet_1.fastq \
#     FASTQ2=tutorial/tutorial_6484_FastqToSam/6484_snippet_2.fastq \
#     OUTPUT=tutorial/output/6484_snippet_fastqtosam.bam \
#     SAMPLE_NAME=NA12878
# 
# # complete example from website
# java -Xmx8G -jar picard.jar FastqToSam \
#     FASTQ=6484_snippet_1.fastq \ #first read file of pair
#     FASTQ2=6484_snippet_2.fastq \ #second read file of pair
#     OUTPUT=6484_snippet_fastqtosam.bam \
#     READ_GROUP_NAME=H0164.2 \ #required; changed from default of A
#     SAMPLE_NAME=NA12878 \ #required
#     LIBRARY_NAME=Solexa-272222 \ #required
#     PLATFORM_UNIT=H0164ALXX140820.2 \
#     PLATFORM=illumina \ #recommended
#     SEQUENCING_CENTER=BI \
#     RUN_DATE=2014-08-20T00:00:00-0400



# --------------
# run FastqToSam
# --------------

# create unaligned BAM files from FASTQ files for scRNA-seq samples

# runtime: up to 1 hour per sample (on laptop)

cd ../FastqToSam

# java -jar ../picard/picard.jar FastqToSam \
#     FASTQ=~/data/gnomex/16030R_single_cell_sequencing/16030R/Fastq/16030X2_HJVMLDMXX/16030X2_HJVMLDMXX_S1_L002_R1_001.fastq.gz \
#     FASTQ2=~/data/gnomex/16030R_single_cell_sequencing/16030R/Fastq/16030X2_HJVMLDMXX/16030X2_HJVMLDMXX_S1_L002_R2_001.fastq.gz \
#     OUTPUT=output/16030X2_HJVMLDMXX_fastqtosam.bam \
#     SAMPLE_NAME=16030X2_HJVMLDMXX

# java -jar ../picard/picard.jar FastqToSam \
#     FASTQ=~/data/gnomex/16030R_single_cell_sequencing/16030R/Fastq/16030X3_HJTWLDMXX/16030X3_HJTWLDMXX_S1_L001_R1_001.fastq.gz \
#     FASTQ2=~/data/gnomex/16030R_single_cell_sequencing/16030R/Fastq/16030X3_HJTWLDMXX/16030X3_HJTWLDMXX_S1_L001_R2_001.fastq.gz \
#     OUTPUT=output/16030X3_HJTWLDMXX_fastqtosam.bam \
#     SAMPLE_NAME=16030X3_HJTWLDMXX

java -jar ../picard/picard.jar FastqToSam \
    FASTQ=~/data/gnomex/16030R_single_cell_sequencing/16030R/Fastq/16030X4_HJTWLDMXX/16030X4_HJTWLDMXX_S2_L001_R1_001.fastq.gz \
    FASTQ2=~/data/gnomex/16030R_single_cell_sequencing/16030R/Fastq/16030X4_HJTWLDMXX/16030X4_HJTWLDMXX_S2_L001_R2_001.fastq.gz \
    OUTPUT=output/16030X4_HJTWLDMXX_fastqtosam.bam \
    SAMPLE_NAME=16030X4_HJTWLDMXX


