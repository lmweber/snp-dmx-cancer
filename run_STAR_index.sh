# Bash script to run STAR to create genome index for STARsolo

# STARsolo is a faster and less buggy alternative to Cell Ranger
# see https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md


# --------------
# generate index
# --------------

# generate STAR genome index using Cell Ranger annotation files (only needs to be run once)

#qsub -V -cwd -pe local 20 -l mem_free=10G,h_vmem=11G,h_fsize=200G run_STAR_index.sh


STAR --runMode genomeGenerate \
--runThreadN 20 \
--genomeDir ../data/STAR_index \
--genomeFastaFiles ../data/GRCh38/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
--sjdbGTFfile ../data/GRCh38/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf \
--genomeSAsparseD 3


