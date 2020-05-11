# Evaluate performance

# evaluate performance of Vireo for demultiplexing samples

# notes:
# ran cellSNP and Vireo using option 1 (not using any information from bulk samples)


#module load conda_R/4.0
#R


library(tidyverse)
library(magrittr)



# -----------------
# load ground truth
# -----------------

# ground truth is stored in sample ID suffixes in cell barcodes

barcodes_all <- read_table("barcodes/barcodes_merged.txt", col_names = "barcode")

df_truth <- barcodes_all %>% 
  mutate(sample_id = gsub("^[A-Za-z]+-", "", barcode))

head(df_truth)



# -----------------
# load Vireo output
# -----------------

# Vireo main outputs are saved in "donor_ids.tsv"

out_vireo <- read_tsv("out_vireo/donor_ids.tsv")

head(out_vireo)

stopifnot(nrow(out_vireo) == nrow(barcodes_all))

# note order of barcodes is different, so need to match them (e.g. using join)
df_truth <- df_truth %>% full_join(out_vireo, by = c("barcode" = "cell"))

df_truth$sample_id %<>% as.factor
df_truth$donor_id %<>% as.factor

head(df_truth)



# ---------------------
# calculate performance
# ---------------------

tbl_truth <- table(truth = df_truth$sample_id, predicted = df_truth$donor_id)
tbl_truth

n_cells <- rowSums(tbl_truth)
n_cells

stopifnot(all(apply(tbl_truth, 1, sum) == n_cells))


# calculate precision and recall

recall <- apply(tbl_truth, 1, function(row) { max(row) / sum(row) })
recall

# note sample order is different
sample_order <- c("donor1", "donor2", "donor0")
precision <- rep(NA, length(sample_order))
for (i in seq_along(sample_order)) {
  precision[i] <- max(tbl_truth[, sample_order[i]] / sum(tbl_truth[, sample_order[i]]))
}
names(precision) <- paste(sample_order, names(recall), sep = "_")
precision


# summary
n_cells
tbl_truth
precision
recall


