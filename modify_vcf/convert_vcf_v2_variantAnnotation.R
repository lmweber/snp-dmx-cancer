# Convert VCF file from genomic coordinates to transcriptomic coordinates

# Steps:
# (1) filter VCF file to keep only variants in 3' UTR regions for faster cellSNP
# runtime (using script "filter_vcf_mod.R")
# (2) manually upload filtered VCF file to Ensembl Variant Effect Predictor
# (VEP) tool to get annotated .vcf / .vep / .txt versions containing both
# genomic and transcriptomic coordinates
# (3) run this script to convert the filtered VCF file from genomic to
# transcriptomic coordinates

# Lukas Weber, 2020-07-18


# This version uses the VariantAnnotation package to create a new VCF file.



# ----------------------
# Compress and index VCF
# ----------------------

# first need to compress and create an index file for the VCF file from VEP, so
# the VariantAnnotation functions can load it

# commands on JHPCE server

# module load htslib
# bgzip genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs_VEP_annotated.vcf
# bcftools index genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs_VEP_annotated.vcf.gz


# --------
# Load VCF
# --------

library(VariantAnnotation)

vcf <- VcfFile(
  "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs_VEP_annotated.vcf.gz", 
  "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs_VEP_annotated.vcf.gz.csi")

vars <- readVcf(
  vcf, 
  param = ScanVcfParam(
    fixed = names(fixed(scanVcfHeader(vcf))), 
    info = names(info(scanVcfHeader(vcf))), 
    geno = names(geno(scanVcfHeader(vcf)))))


# ---------------------------------------
# Load transcriptomic annotation from VEP
# ---------------------------------------

# load annotated TXT file from VEP containing transcript ID and position
file_vep_txt <- "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs_VEP_annotated.txt"
vep_txt <- read.delim(file_vep_txt)

# get transcript ID and position
transcript_id <- vep_txt$Feature
transcript_pos <- vep_txt$cDNA_position


# ----------------------------------------------------
# Create new VCF object containing transcriptomic info
# ----------------------------------------------------

# create new rowRanges object for VCF
rowranges <- GRanges(
  Rle(transcript_id), 
  ranges(rowRanges(vars)), 
  strand(rowRanges(vars)), 
  mcols(rowRanges(vars)))

rowranges

# filter out variants with missing transcript position
ix_filter <- transcript_pos == "-"
table(ix_filter)
mean(ix_filter)  ## filters out 17.5% of variants

rowranges <- rowranges[!ix_filter, ]
rowranges

# replace position with transcriptomic position
transcript_pos_filt <- transcript_pos[!ix_filter]
rr <- ranges(rowranges)
start(rr) <- as.integer(transcript_pos_filt)
end(rr) <- as.integer(transcript_pos_filt)

ranges(rowranges) <- rr
rowranges

# create VCF object with new rowRanges
vars <- vars[!ix_filter, ]
vars

rowRanges(vars) <- rowranges

rowRanges(vars)


# -----------------
# Save new VCF file
# -----------------

writeVcf(vars, "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs_transcriptomic_v2.vcf.gz", index = TRUE)

