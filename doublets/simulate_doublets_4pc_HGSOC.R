# -----------------
# Simulate doublets
# -----------------

# simulate doublets by combining some percentage of cell barcodes

# notes:
# - generate a lookup table to replace and combine some percentage of cell barcodes
# - parse through merged BAM file and merged barcodes file using sed to replace cell barcodes
# - then continue with cellSNP and Vireo and evaluate demultiplexing performance

# runtime: 3 hours


# interactive session on cluster
#qrsh -l mem_free=10G,h_vmem=20G,h_fsize=100G


# ----------
# Parameters
# ----------

# proportion of doublets to simulate (i.e. proportion of final number of cells)
prop_doublets <- 0.04
# number of cells to replace (taking into account reduction in final number of cells)
prop_doublets_corrected <- prop_doublets * (1 - prop_doublets)


# --------------------------
# Generate replacement table
# --------------------------

library(tidyverse)
library(magrittr)

# load merged list of cell barcodes
df_barcodes <- read_table("../outputs/cellranger_testing_HGSOC/barcodes_merged.tsv", 
                          col_names = "barcode")

# number of cells
n_cells <- nrow(df_barcodes)
# number of cells to combine with other cells as doublets
n_doublets <- round(prop_doublets_corrected * n_cells)

n_cells
n_doublets

# select random sets of (i) cells to replace and (ii) replacements
set.seed(123)
ix_original <- sample(seq_len(n_cells), n_doublets)
ix_replacement <- sample(setdiff(seq_len(n_cells), ix_original), n_doublets)

# generate lookup table / replacement table for cell barcodes
df_lookup <- tibble(
  original = df_barcodes$barcode[ix_original], 
  replacement = df_barcodes$barcode[ix_replacement]
)

# save lookup table
filename <- paste0("replacement_table_doublets_", prop_doublets * 100, "pc_HGSOC.tsv")
write_tsv(df_lookup, paste0("../outputs/doublets/", filename))


# ------------------------------
# Generate updated barcodes file
# ------------------------------

# generate updated list of barcodes
barcodes_merged_new <- df_barcodes$barcode
barcodes_merged_new[ix_original] <- df_barcodes$barcode[ix_replacement]
stopifnot(length(unique(barcodes_merged_new)) == n_cells - n_doublets)
stopifnot(length(barcodes_merged_new) == n_cells)

# save barcodes file
filename <- paste0("../outputs/doublets/barcodes_merged_", prop_doublets * 100, "pc_HGSOC.tsv")
write_tsv(tibble(barcodes_merged_new), filename, col_names = FALSE)


# --------------------
# Generate sed command
# --------------------

# convert lookup table into a multiple replacement sed command for parsing through BAM file

cmd1 <- paste0("s/", df_lookup$original, "/", df_lookup$replacement, "/g")
cmd2 <- paste0(cmd1, collapse = ";")
cmd_bam <- paste0("sed -e '", cmd2, "'")
#cmd_bam <- paste0(cmd3, " ../outputs/HGSOC/BAM_merged/BAM_merged.bam")

str(cmd_bam)


# combine with "samtools view" commands and pipes to convert BAM to SAM and back again
# see script "parse_BAM_files.sh" in main pipeline for more details

cmd_full <- paste0(
  "samtools view -h ../outputs/HGSOC/BAM_merged/BAM_merged.bam | ", 
  cmd_bam, " | ", 
  "samtools view -bo ../outputs/HGSOC/doublets/BAM_merged_doublets_", prop_doublets * 100, "pc_HGSOC.bam"
)

str(cmd_full)
cmd_full

# save a shell script
writeLines(cmd_full, "doublets/doublets_BAM_4pc_HGSOC.sh")

