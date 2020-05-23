# Evaluate performance

# evaluate performance of Vireo for demultiplexing samples

# notes:
# ran cellSNP and Vireo using option 1 (not using any information from bulk samples)


#module load conda_R/4.0
#R


library(tidyverse)
library(magrittr)
library(Matrix)



# -----------------
# load ground truth
# -----------------

# ground truth is stored in sample ID suffixes in cell barcodes

barcodes_all <- read_table("../CellRanger/barcodes/barcodes_merged.txt", col_names = "barcode")

df_truth <- barcodes_all %>% 
  mutate(sample_id = gsub("^[A-Za-z]+-", "", barcode))

head(df_truth)



# -----------------
# load Vireo output
# -----------------

# Vireo main outputs are saved in "donor_ids.tsv"

out_vireo <- read_tsv("../CellRanger/out_vireo/donor_ids.tsv")

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



# ----------------------------------------
# investigate unassigned and doublet cells
# ----------------------------------------

# note Cell Ranger output files have been gunzipped

# details on how to load Cell Ranger output matrices into R:
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices


df_truth <- df_truth %>% 
  mutate(unassigned = donor_id == "unassigned") %>% 
  mutate(doublet = donor_id == "doublet")

df_truth


# load Cell Ranger output

dir_out_X2 <- "../CellRanger/16030X2_HJVMLDMXX/outs/filtered_feature_bc_matrix"
dir_out_X3 <- "../CellRanger/16030X3_HJTWLDMXX/outs/filtered_feature_bc_matrix"
dir_out_X4 <- "../CellRanger/16030X4_HJTWLDMXX/outs/filtered_feature_bc_matrix"


# code from:
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices

mat_X2 <- readMM(file.path(dir_out_X2, "matrix.mtx.gz"))
feature_names_X2 <- read.delim(file.path(dir_out_X2, "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
barcode_names_X2 <- read.delim(file.path(dir_out_X2, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
colnames(mat_X2) <- gsub("-1$", "-X2", barcode_names_X2$V1)
rownames(mat_X2) <- feature_names_X2$V1

mat_X3 <- readMM(file.path(dir_out_X3, "matrix.mtx.gz"))
feature_names_X3 <- read.delim(file.path(dir_out_X3, "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
barcode_names_X3 <- read.delim(file.path(dir_out_X3, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
colnames(mat_X3) <- gsub("-1$", "-X3", barcode_names_X3$V1)
rownames(mat_X3) <- feature_names_X3$V1

mat_X4 <- readMM(file.path(dir_out_X4, "matrix.mtx.gz"))
feature_names_X4 <- read.delim(file.path(dir_out_X4, "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
barcode_names_X4 <- read.delim(file.path(dir_out_X4, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
colnames(mat_X4) <- gsub("-1$", "-X4", barcode_names_X4$V1)
rownames(mat_X4) <- feature_names_X4$V1


# check
mat_X2[1:6, 1:6]


# calculate number of genes and proportion mitochondrial reads per cell

is_mito_X2 <- grepl("^MT-", feature_names_X2$V2)
is_mito_X3 <- grepl("^MT-", feature_names_X3$V2)
is_mito_X4 <- grepl("^MT-", feature_names_X4$V2)

table(is_mito_X2)
table(is_mito_X3)
table(is_mito_X4)

prop_mito_X2 <- colSums(mat_X2[is_mito_X2, ]) / colSums(mat_X2)
prop_mito_X3 <- colSums(mat_X3[is_mito_X3, ]) / colSums(mat_X3)
prop_mito_X4 <- colSums(mat_X4[is_mito_X4, ]) / colSums(mat_X4)

n_genes_X2 <- colSums(as.matrix(mat_X2) > 0)
n_genes_X3 <- colSums(as.matrix(mat_X3) > 0)
n_genes_X4 <- colSums(as.matrix(mat_X4) > 0)


# combine and check order
n_genes_all <- c(n_genes_X2, n_genes_X3, n_genes_X4)
stopifnot(length(n_genes_all) == nrow(df_truth))
n_genes_all <- n_genes_all[df_truth$barcode]
stopifnot(all(names(n_genes_all) == df_truth$barcode))

prop_mito_all <- c(prop_mito_X2, prop_mito_X3, prop_mito_X4)
stopifnot(length(prop_mito_all) == nrow(df_truth))
prop_mito_all <- prop_mito_all[df_truth$barcode]
stopifnot(all(names(prop_mito_all) == df_truth$barcode))

# add to table
df_truth$n_genes <- n_genes_all
df_truth$prop_mito <- prop_mito_all


# results: number of genes and proportion mitochondrial reads per cell,
# stratified by unassigned and doublet status from Vireo

df_truth %>% 
  group_by(sample_id, unassigned) %>% 
  summarize(n_cells = n(), 
            mean_n_genes = mean(n_genes), 
            mean_prop_mito = mean(prop_mito))

df_truth %>% 
  group_by(sample_id, doublet) %>% 
  summarize(n_cells = n(), 
            mean_n_genes = mean(n_genes), 
            mean_prop_mito = mean(prop_mito))



