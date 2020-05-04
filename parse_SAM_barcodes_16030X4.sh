# Parse SAM files from Cell Ranger to add unique sample IDs to cell barcodes

# notes:
# syntax to search and replace: sed -i "s/regexp/replacement/g"
# regular expression matches cell barcode with format "CB:Z:AGCTTCCAGCAGTCTT-1"
# where the "-1" is the default sample ID added by Cell Ranger
# then we replace "-1" with a unique sample ID

# see also https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam


#qsub -V -cwd -l mem_free=100G,h_vmem=200G,h_fsize=300G parse_SAM_barcodes.sh


sed -i "s/\(CB\:Z\:[A-Z]\+\)\-1/\1\-X2/g" 16030X4_HJTWLDMXX/outs/possorted_genome_bam.sam


