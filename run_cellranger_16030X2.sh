# Bash script to run cellranger to create BAM files for scRNA-seq samples

# see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct
# and https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/job-submission-mode


########### note: does not work on our JHPCE cluster, possibly due to memory allocation issues on cluster
########### switch to using STARsolo instead (see other scripts)
########### this script might work on a standalone linux system though


#qsub -V -cwd -pe local 10 -l mem_free=5G,h_vmem=10G run_cellranger.sh


cellranger count --id=16030X2_HJVMLDMXX \
--fastqs=../data/16030R/Fastq/16030X2_HJVMLDMXX \
--sample=16030X2_HJVMLDMXX \
--transcriptome=../data/GRCh38/refdata-cellranger-GRCh38-3.0.0 \
--expect-cells=5000 \
--jobmode=local \
--localcores=10 \
--localmem=45

# cellranger count --id=16030X3_HJTWLDMXX \
# --fastqs=../data/16030R/Fastq/16030X3_HJTWLDMXX \
# --sample=16030X3_HJTWLDMXX \
# --transcriptome=../data/GRCh38/refdata-cellranger-GRCh38-3.0.0 \
# --expect-cells=5000 \
# --jobmode=local \
# --localcores=10 \
# --localmem=45

# cellranger count --id=16030X4_HJTWLDMXX \
# --fastqs=../data/16030R/Fastq/16030X4_HJTWLDMXX \
# --sample=16030X4_HJTWLDMXX \
# --transcriptome=../data/GRCh38/refdata-cellranger-GRCh38-3.0.0 \
# --expect-cells=5000 \
# --jobmode=local \
# --localcores=10 \
# --localmem=45


