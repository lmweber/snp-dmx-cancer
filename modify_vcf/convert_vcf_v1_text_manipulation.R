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


# The following code loads the VCF file as text, breaks it up into header /
# column names / data, modifies the CHROM and POS data columns to contain
# transcript ID and position instead of chromosome ID and position, and outputs
# as a new VCF file. (It may also be possible to do this using the
# VariantAnnotation package, but simple text manipulation also seems to work.)


# --------
# Load VCF
# --------

# load annotated VCF file as text
vcf_file <- "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs_VEP_annotated.vcf"
vcf_lines <- readLines(vcf_file)

# extract header (lines beginning with ##)
ix_header <- grep("^##", vcf_lines)
ix_header_logi <- grepl("^##", vcf_lines)
vcf_header <- vcf_lines[ix_header_logi]

# extract column names in VCF format
vcf_colnames <- vcf_lines[max(ix_header) + 1]

# load main part of VCF as data frame
library(readr)
vcf_df <- read.delim(textConnection(vcf_lines[!ix_header_logi]))

head(vcf_df)


# ---------------------------------------------------------
# Load annotated TXT output file from VEP and match entries
# ---------------------------------------------------------

# load annotated TXT file from VEP containing transcript ID and position
file_vep_txt <- "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs_VEP_annotated.txt"
vep_txt <- read.delim(file_vep_txt)

# check rows match
nrow(vcf_df) == nrow(vep_txt)
all(vcf_df$ID == vep_txt$X.Uploaded_variation)

# replace chromosome ID and position with transcript ID and position in VCF file
transcript_id <- vep_txt$Feature
transcript_pos <- vep_txt$cDNA_position

vcf_df_mod <- vcf_df

vcf_df_mod$X.CHROM <- transcript_id
vcf_df_mod$POS <- transcript_pos

head(vcf_df_mod)

# filter out variants with missing transcript position
ix_filter <- transcript_pos == "-"
table(ix_filter)
mean(ix_filter)  ## filters out 17.5% of variants

vcf_df_mod <- vcf_df_mod[!ix_filter, ]

dim(vcf_df_mod)
head(vcf_df_mod)


# ------------------------
# Output modified VCF file
# ------------------------

# save modified VCF by putting together parts (header, column names, data)
vcf_file_out <- "../../data/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.threeUTRs_transcriptomic.vcf"

writeLines(c(vcf_header, vcf_colnames), vcf_file_out)

write.table(vcf_df_mod, file = vcf_file_out, append = TRUE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

