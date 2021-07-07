###############################################################
# Script to subset VCF files using list of SNPs from MEGA array
# Lukas Weber, June 2021
###############################################################

# subset SNPs from MEGA array in VCF files from bcftools (bulk RNA-seq samples)
# and 1000 Genomes; SNPs are identified by chromosome and position

library(readr)


# ----------
# load files
# ----------

# read in VCF file from bcftools (bulk RNA-seq samples)

vcf_bcftools <- read_tsv("../../genotype/bcftools/bcftools_HGSOC_rehead.vcf", skip = 221)
head(vcf_bcftools)
nrow(vcf_bcftools)  # 605,375

# read in VCF file from 1000 Genomes Project (provided by cellSNP authors)

vcf_1000GenomesUnfilt <- read_tsv("../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf", skip = 234)
head(vcf_1000GenomesUnfilt)
nrow(vcf_1000GenomesUnfilt)  # 7,416,067

# read in VCF file from 1000 Genomes Project (filtered to 3' UTR region)

vcf_1000GenomesFilt <- read_tsv("../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf", skip = 225)
head(vcf_1000GenomesFilt)
nrow(vcf_1000GenomesFilt)  # 85,594

# read in BED file from MEGA SNP array (Infinium Multi-Ethnic Global-8 v1.0) (build 38)
# https://support.illumina.com/array/array_kits/infinium-multi-ethnic-global-8-kit/downloads.html

bed_MEGA <- read_tsv("../../data/MEGA/Multi-EthnicGlobal_D2.bed", skip = 1, col_names = FALSE)
head(bed_MEGA)
nrow(bed_MEGA)  # 1,733,365


# ----------
# match SNPs
# ----------

# format into consistent strings

snps_bcftools <- paste0(vcf_bcftools[, 1, drop = TRUE], "_", vcf_bcftools[, 2, drop = TRUE])
length(snps_bcftools)

snps_1000GenomesUnfilt <- paste0("chr", vcf_1000GenomesUnfilt[, 1, drop = TRUE], "_", vcf_1000GenomesUnfilt[, 2, drop = TRUE])
length(snps_1000GenomesUnfilt)

snps_1000GenomesFilt <- paste0("chr", vcf_1000GenomesFilt[, 1, drop = TRUE], "_", vcf_1000GenomesFilt[, 2, drop = TRUE])
length(snps_1000GenomesFilt)

snps_MEGA <- paste0(bed_MEGA[, 1, drop = TRUE], "_", bed_MEGA[, 3, drop = TRUE])
length(snps_MEGA)


# calculate overlapping set sizes

head(snps_bcftools)
head(snps_1000GenomesUnfilt)
head(snps_1000GenomesFilt)
head(snps_MEGA)

sum(snps_MEGA %in% snps_bcftools)  # 45846 out of 605375 (7.6%)
sum(snps_MEGA %in% snps_1000GenomesUnfilt)  # 638665 out of 7416067 (8.6%)
sum(snps_MEGA %in% snps_1000GenomesFilt)  # 14021 out of 85594 (16.4%)

sum(snps_bcftools %in% snps_MEGA)  # 45846 out of 1733365 (2.6%)
sum(snps_1000GenomesUnfilt %in% snps_MEGA)  # 638668 out of 1733365 (3.7%)
sum(snps_1000GenomesFilt %in% snps_MEGA)  # 14150 out of 1733365 (0.8%)


# also calculate previous overlap

sum(snps_bcftools %in% snps_1000GenomesUnfilt)  # 356814 out of 7416067 (4.8%)
sum(snps_bcftools %in% snps_1000GenomesFilt)  # 28502 out of 85594 (33.3%)


# --------------
# save VCF files
# --------------

# save new VCF files containing only the overlapping SNPs

ix_bcftools <- sample(which(snps_bcftools %in% snps_MEGA))  # note: permute to avoid errors when running cellSNP
ix_1000GenomesUnfilt <- which(snps_1000GenomesUnfilt %in% snps_MEGA)
ix_1000GenomesFilt <- which(snps_1000GenomesFilt %in% snps_MEGA)

file_bcftools <- "../../genotype/MEGA/subset_MEGA_bcftools.vcf"
file_1000GenomesUnfilt <- "../../genotype/MEGA/subset_MEGA_1000GenomesUnfilt.vcf"
file_1000GenomesFilt <- "../../genotype/MEGA/subset_MEGA_1000GenomesFilt.vcf"

system(paste0("head -n 222 ../../genotype/bcftools/bcftools_HGSOC_rehead.vcf > ", file_bcftools))
write_tsv(vcf_bcftools[ix_bcftools, ], file = file_bcftools, append = TRUE)

system(paste0("head -n 235 ../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf > ", file_1000GenomesUnfilt))
write_tsv(vcf_1000GenomesUnfilt[ix_1000GenomesUnfilt, ], file = file_1000GenomesUnfilt, append = TRUE)

system(paste0("head -n 226 ../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf > ", file_1000GenomesFilt))
write_tsv(vcf_1000GenomesFilt[ix_1000GenomesFilt, ], file = file_1000GenomesFilt, append = TRUE)

