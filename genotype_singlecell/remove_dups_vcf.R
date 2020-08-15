###################################################
# R script to remove duplicate variants in VCF file
###################################################

# R script to remove duplicate variants in merged VCF file from the previous
# step (concatenate_vcf.sh). Duplicate variants are those that are found in
# multiple samples, so removing them is likely to improve sample demultiplexing
# performance in subsequent steps.


# module load conda_R/4.0
# Rscript remove_dups_vcf.R


# --------
# load VCF
# --------

# load annotated VCF file as text
vcf_file <- "../../singlecell/cellSNP_singlecell_merged/cellSNP.cells-merged.vcf"
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


# --------------------------
# identify duplicate entries
# --------------------------

# duplicate entries are those with the same chromosome (column 1) and position (column 2)

entries <- paste(vcf_df[, 1], vcf_df[, 2], sep = "_")

ix_dups <- duplicated(entries)
ix_dups_all <- entries %in% entries[ix_dups]

table(ix_dups)
table(ix_dups_all)

vcf_df_nodups <- vcf_df[!ix_dups_all, , drop = FALSE]


# -------------
# save VCF file
# -------------

# save VCF file without duplicate entries

vcf_file_out <- "../../singlecell/cellSNP_singlecell_merged/cellSNP.cells-merged-nodups.vcf"

writeLines(c(vcf_header, vcf_colnames), vcf_file_out)

write.table(vcf_df_nodups, file = vcf_file_out, append = TRUE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# ---------------------
# save gzipped VCF file
# ---------------------

# also save VCF file in gzipped format

system("gzip -c ../../singlecell/cellSNP_singlecell_merged/cellSNP.cells-merged-nodups.vcf > ../../singlecell/cellSNP_singlecell_merged/cellSNP.cells-merged-nodups.vcf.gz")

