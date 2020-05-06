# Run cellSNP to genotype cells

# see https://vireosnp.readthedocs.io/en/latest/genotype.html
# and https://github.com/single-cell-genetics/cellSNP


# notes:
# requires merged BAM file and merged barcodes file, each containing unique sample IDs
# run cellSNP in mode 1 (using .vcf file containing common variants)


#qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=30G,h_fsize=300G run_cellSNP.sh


cellSNP \
-s bam_merged/bam_merged.bam \
-b barcodes/barcodes_merged.txt \
-O out_cellSNP \
-R ../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
-p 10 \
--minMAF 0.1 \
--minCOUNT 20


