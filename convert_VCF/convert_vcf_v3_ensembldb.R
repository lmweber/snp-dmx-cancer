# Convert VCF file from genomic coordinates to transcriptomic coordinates

# Steps:
# (1) filter VCF file to keep only variants in 3' UTR regions for faster cellSNP
# runtime (using script "filter_vcf_mod.R")
# (2) convert genomic coordinates to transcriptomic coordinates in VCF file
# using ensembldb package

# Lukas Weber, 2020-07-19


# This version uses the ensembldb package to convert genomic to transcriptomic
# coordinates.



# -------------
# Load VCF file
# -------------

# note header info in vcf file refers to b37/hg19 genome build, but the data in
# the vcf file has been lifted over to hg38 build
# see https://sourceforge.net/projects/cellsnp/files/SNPlist/

# this requires some hacking of the vcf object further below

library(VariantAnnotation)

vcf <- VcfFile(
  "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf.bgz", 
  "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf.bgz.tbi")

vars <- readVcf(
  vcf, 
  param = ScanVcfParam(
    fixed = names(fixed(scanVcfHeader(vcf))), 
    info = names(info(scanVcfHeader(vcf))), 
    geno = names(geno(scanVcfHeader(vcf)))))


# -------------------
# Load Ensembldb info
# -------------------

library(AnnotationHub)
library(ensembldb)

ah <- AnnotationHub()
EnsDb.Hsapiens.v94 <- ah[["AH64923"]]
edb <- EnsDb.Hsapiens.v94

# get seqinfo
seqinfo_new <- seqinfo(edb)


# ---------------------
# Manipulate vcf object
# ---------------------

# manipulate vcf object to update build info (see comments above); keep only
# levels for chromosomes 1-22 and X

# this is not an ideal way to do this, but seems to work

str(vars)

vars@rowRanges@seqnames@values <- droplevels(vars@rowRanges@seqnames@values)
vars@rowRanges@seqinfo <- seqinfo_new[c(1:22, "X"), ]


# --------------------------------------
# Convert genomic to transcriptomic info
# --------------------------------------

rr <- rowRanges(vars)

rr_tx <- genomeToTranscript(rr, edb)

rr_new <- unlist(rr_tx)


# then put new rowRanges back into VCF object

# however note there are multiple transcript IDs for some variants
# not clear how to write this back into VCF file

