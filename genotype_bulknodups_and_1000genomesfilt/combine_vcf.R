###############################
# R script to combine VCF files
###############################

# combine variants from:
# (i) 1000 Genomes Project (filtered)
# (ii) bulk genotyping (no duplicates)


# module load conda_R/4.0
# Rscript remove_dups_vcf.R


# -------------------------------------------
# load VCF 1: 1000 Genomes Project (filtered)
# -------------------------------------------

# load annotated VCF file as text
vcf_file_1 <- "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs.vcf"
vcf_lines_1 <- readLines(vcf_file_1)

# extract header (lines beginning with ##)
ix_header_1 <- grep("^##", vcf_lines_1)
ix_header_logi_1 <- grepl("^##", vcf_lines_1)
vcf_header_1 <- vcf_lines_1[ix_header_logi_1]

# extract column names in VCF format
vcf_colnames_1 <- vcf_lines_1[max(ix_header_1) + 1]

# load main part of VCF as data frame
library(readr)
vcf_df_1 <- read.delim(textConnection(vcf_lines_1[!ix_header_logi_1]))

head(vcf_df_1)


# ---------------------------------------------------
# load VCF 2: bulk samples genotyping (no duplicates)
# ---------------------------------------------------

# load annotated VCF file as text
vcf_file_2 <- "../../genotype_bulk/cellSNP_bulk_merged/cellSNP.cells-merged-nodups.vcf"
vcf_lines_2 <- readLines(vcf_file_2)

# extract header (lines beginning with ##)
ix_header_2 <- grep("^##", vcf_lines_2)
ix_header_logi_2 <- grepl("^##", vcf_lines_2)
vcf_header_2 <- vcf_lines_2[ix_header_logi_2]

# extract column names in VCF format
vcf_colnames_2 <- vcf_lines_2[max(ix_header_2) + 1]

# load main part of VCF as data frame
library(readr)
vcf_df_2 <- read.delim(textConnection(vcf_lines_2[!ix_header_logi_2]))

head(vcf_df_2)


# ---------------
# merge VCF files
# ---------------

# convert column 1 to same format
vcf_df_2[, 1] <- gsub("^chr", "", vcf_df_2[, 1])

# merge
cols_keep <- 1:8
vcf_df_merged <- rbind(vcf_df_1, vcf_df_2[, cols_keep])

dim(vcf_df_1)
dim(vcf_df_2)
dim(vcf_df_merged)


# -------------
# save VCF file
# -------------

vcf_file_out <- "../../genotype_bulknodups_and_1000genomesfilt/cellSNP.cells-combined.vcf"

writeLines(c(vcf_header_1, vcf_colnames_1), vcf_file_out)

write.table(vcf_df_merged, file = vcf_file_out, append = TRUE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# ---------------------
# save gzipped VCF file
# ---------------------

# also save in gzipped format

system("gzip -c ../../genotype_bulknodups_and_1000genomesfilt/cellSNP.cells-combined.vcf > ../../genotype_bulknodups_and_1000genomesfilt/cellSNP.cells-combined.vcf.gz")

