# Construct a VCF file for use cellSNP by intersecting the cellSNP-recommended
# VCF file with 3' UTR regions.

# Peter Hickey
# 2020-03-04

# modified by Lukas Weber for use on JHU JHPCE cluster
# 2020-07-09


# load R module
#module load conda_R


# Sort and index the cellSNP-recommended VCF file ------------------------------

# note: bcftools installed locally and available via PATH

cmd <- paste0(
  "gunzip ", 
  "-c ../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz ", 
  "| ", 
  "bcftools sort ", 
  "-o ../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.bgz -O z")
system(cmd)

cmd <- "bcftools index ../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.bgz"
system(cmd)


# Construct 3' UTRs ------------------------------------------------------------

library(AnnotationHub)
library(ensembldb)

ah <- AnnotationHub()
EnsDb.Hsapiens.v94 <- ah[["AH64923"]]
three_utrs <- threeUTRsByTranscript(EnsDb.Hsapiens.v94)
three_utrs <- keepSeqlevels(three_utrs, c(1:22, "X"), pruning.mode = "coarse")
reduced_three_utrs <- reduce(unlist(three_utrs))


# Load VCF for loci in 3' UTRs -------------------------------------------------

library(VariantAnnotation)

vcf <- VcfFile(
  "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.bgz",
  "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.bgz.csi")

vars <- readVcf(
  vcf,
  param = ScanVcfParam(
    fixed = names(fixed(scanVcfHeader(vcf))),
    info = names(info(scanVcfHeader(vcf))),
    geno = names(geno(scanVcfHeader(vcf))),
    which = reduced_three_utrs))


# Write to disk as VCF file ----------------------------------------------------

writeVcf(
  vars, 
  "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf", 
  index = TRUE)


# Unzip VCF file ---------------------------------------------------------------

cmd <- paste0(
  "gunzip -c ", 
  "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf.bgz ", 
  "> ", 
  "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf")
system(cmd)

