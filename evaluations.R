# Evaluate performance

# evaluate performance of Vireo for demultiplexing samples

# notes:
# ran cellSNP and Vireo using option 1 (not using any information from bulk samples)


#module load conda_R/4.0
#R


library(tidyverse)
library(magrittr)
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(scran)



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



# ---------------------
# precision-recall plot
# ---------------------

df_pr <- DataFrame(
  precision = precision, 
  recall = recall, 
  sample = factor(names(recall))
)

p <- 
  as.data.frame(df_pr) %>% 
  ggplot(aes(x = recall, y = precision, color = sample)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = unname(palette.colors(palette = "Okabe-Ito"))) + 
  xlim(c(0.95, 1)) + 
  ylim(c(0.95, 1)) + 
  coord_fixed() + 
  theme_bw() + 
  theme(panel.grid = element_blank())

p

ggsave("../plots/precision_recall.pdf", width = 4, height = 3.25)



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



# --------------------
# calculate QC metrics
# --------------------

# create SCE object

sample_ids_all <- factor(c(rep("X2", ncol(mat_X2)), rep("X3", ncol(mat_X3)), rep("X4", ncol(mat_X4))))
barcodes_all <- c(colnames(mat_X2), colnames(mat_X3), colnames(mat_X4))

n_cells <- sum(ncol(mat_X2), ncol(mat_X3), ncol(mat_X4))
stopifnot(length(sample_ids_all) == n_cells)
stopifnot(length(barcodes_all) == n_cells)

col_data <- DataFrame(
  sample_id = sample_ids_all, 
  barcode = barcodes_all
)

stopifnot(all(rownames(mat_X2) == rownames(mat_X3)))
stopifnot(all(rownames(mat_X2) == rownames(mat_X4)))
stopifnot(nrow(mat_X2) == nrow(mat_X3))
stopifnot(nrow(mat_X2) == nrow(mat_X4))

row_data <- DataFrame(
  gene_id = feature_names_X2$V1, 
  symbol = feature_names_X2$V2
)

counts_all <- cbind(mat_X2, mat_X3, mat_X4)

sce <- SingleCellExperiment(
  assays = list(counts = counts_all), 
  rowData = row_data, 
  colData = col_data
)

sce
rowData(sce)
colData(sce)


# calculate QC metrics using Bioconductor functions

is_mito <- grepl("^MT-", rowData(sce)$symbol)
table(is_mito)

rowData(sce)$is_mito <- is_mito

# calculate per cell QC metrics
sce <- addPerCellQC(sce, subsets = list(mito = is_mito))

colData(sce)

# identify outliers (3 MADs, using log scale for "lower")
qc_lib <- isOutlier(colData(sce)$sum, log = TRUE, type = "lower")
qc_nexprs <- isOutlier(colData(sce)$detected, log = TRUE, type = "lower")
qc_mito <- isOutlier(colData(sce)$subsets_mito_percent, type = "higher")

# note: thresholds are not appropriate due to large number of dead cells
attr(qc_lib, "thresholds")
attr(qc_nexprs, "thresholds")
attr(qc_mito, "thresholds")


# alternatively: calculate metrics separately for each sample (one batch per sample)
qc_lib2 <- isOutlier(colData(sce)$sum, log = TRUE, type = "lower", batch = colData(sce)$sample_id)
qc_nexprs2 <- isOutlier(colData(sce)$detected, log = TRUE, type = "lower", batch = colData(sce)$sample_id)
qc_mito2 <- isOutlier(colData(sce)$subsets_mito_percent, type = "higher", batch = colData(sce)$sample_id)

# thresholds differ substantially between samples, but are still not appropriate
attr(qc_lib2, "thresholds")
attr(qc_nexprs2, "thresholds")
attr(qc_mito2, "thresholds")


# alternatively: specify simple thresholds instead
quantile(colData(sce)$sum, seq(0, 1, by = 0.1))
quantile(colData(sce)$detected, seq(0, 1, by = 0.1))
quantile(colData(sce)$subsets_mito_percent, seq(0, 1, by = 0.1))

# simple thresholds chosen based on distributions above
thresh_lib <- 2000
thresh_detected <- 500
thresh_mito_pc <- 30

# manual QC using thresholds
qc_lib_thresh <- colData(sce)$sum < thresh_lib
qc_nexprs_thresh <- colData(sce)$detected < thresh_detected
qc_mito_thresh <- colData(sce)$subsets_mito_percent < thresh_mito_pc

discard <- qc_lib_thresh | qc_nexprs_thresh | qc_mito_thresh

DataFrame(
  LibSize = sum(qc_lib_thresh), 
  NExprs = sum(qc_nexprs_thresh), 
  MitoProp = sum(qc_mito_thresh), 
  TotalFiltered = sum(discard), 
  NCells = length(discard)
)



# --------
# QC plots
# --------

# combine data frames

stopifnot(length(discard) == ncol(sce))
colData(sce) <- cbind(colData(sce), qc_lib_thresh, qc_nexprs_thresh, qc_mito_thresh, discard)

stopifnot(nrow(df_truth) == ncol(sce))
stopifnot(all(df_truth$barcode == colData(sce)$barcode))
stopifnot(all(df_truth$sample_id == colData(sce)$sample_id))
colData(sce) <- cbind(colData(sce), df_truth[, -c(1:2)])

colData(sce)


# plot QC thresholds

plotColData(sce, x = "detected", y = "subsets_mito_percent", colour_by = "sample_id")

plotColData(sce, x = "detected", y = "subsets_mito_percent", colour_by = "qc_lib_thresh")
plotColData(sce, x = "detected", y = "subsets_mito_percent", colour_by = "qc_nexprs_thresh")
plotColData(sce, x = "detected", y = "subsets_mito_percent", colour_by = "qc_mito_thresh")

plotColData(sce, x = "detected", y = "subsets_mito_percent", colour_by = "discard")

table(colData(sce)$discard)
mean(colData(sce)$discard)  ## discarding very large proportion of cells


# calculate PCA (with/without QC filtered cells)

sce_sub <- sce[, !discard]
dim(sce_sub)


# calculate PCA and UMAP: with QC filtered cells still included
set.seed(100)
clus <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster = clus)
sce <- logNormCounts(sce)
sce

dec <- modelGeneVar(sce)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
str(top_hvgs)

set.seed(100)
sce <- runPCA(sce, subset_row = top_hvgs)
reducedDimNames(sce)
dim(reducedDim(sce, "PCA"))

set.seed(100)
sce <- runUMAP(sce, dimred = "PCA")


# calculate PCA: without QC filtered cells
set.seed(100)
clus <- quickCluster(sce_sub)
sce_sub <- computeSumFactors(sce_sub, cluster = clus)
sce_sub <- logNormCounts(sce_sub)
sce_sub

dec <- modelGeneVar(sce_sub)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
str(top_hvgs)

set.seed(100)
sce_sub <- runPCA(sce_sub, subset_row = top_hvgs)
reducedDimNames(sce_sub)
dim(reducedDim(sce_sub, "PCA"))

set.seed(100)
sce_sub <- runUMAP(sce_sub, dimred = "PCA")



# ------------------------------------------------
# PCA plots: with QC filtered cells still included
# ------------------------------------------------

plotReducedDim(sce, dimred = "PCA", colour_by = "sample_id")

plotReducedDim(sce, dimred = "PCA", colour_by = "qc_lib_thresh")
plotReducedDim(sce, dimred = "PCA", colour_by = "qc_nexprs_thresh")
plotReducedDim(sce, dimred = "PCA", colour_by = "qc_mito_thresh")

plotReducedDim(sce, dimred = "PCA", colour_by = "discard")

# plot Vireo unassigned/doublets in PCA space
plotReducedDim(sce, dimred = "PCA", colour_by = "unassigned")
plotReducedDim(sce, dimred = "PCA", colour_by = "doublet")


# UMAP
plotReducedDim(sce, dimred = "UMAP", colour_by = "sample_id")
plotReducedDim(sce, dimred = "UMAP", colour_by = "qc_lib_thresh")
plotReducedDim(sce, dimred = "UMAP", colour_by = "qc_nexprs_thresh")
plotReducedDim(sce, dimred = "UMAP", colour_by = "qc_mito_thresh")
plotReducedDim(sce, dimred = "UMAP", colour_by = "discard")
plotReducedDim(sce, dimred = "UMAP", colour_by = "n_genes")
plotReducedDim(sce, dimred = "UMAP", colour_by = "prop_mito")
plotReducedDim(sce, dimred = "UMAP", colour_by = "unassigned")
plotReducedDim(sce, dimred = "UMAP", colour_by = "doublet")


# plot QC metrics for Vireo unassigned/doublets

# unassigned
plotColData(sce, x = "unassigned", y = "sum", colour_by = "qc_lib_thresh") + scale_y_log10()
plotColData(sce, x = "unassigned", y = "detected", colour_by = "qc_nexprs_thresh") + scale_y_log10()
plotColData(sce, x = "unassigned", y = "subsets_mito_percent", colour_by = "qc_mito_thresh")

# doublets
plotColData(sce, x = "doublet", y = "sum", colour_by = "qc_lib_thresh") + scale_y_log10()
plotColData(sce, x = "doublet", y = "detected", colour_by = "qc_nexprs_thresh") + scale_y_log10()
plotColData(sce, x = "doublet", y = "subsets_mito_percent", colour_by = "qc_mito_thresh")



# ------------------------------------
# PCA plots: without QC filtered cells
# ------------------------------------

# first PC corresponds to sample ID
plotReducedDim(sce_sub, dimred = "PCA", colour_by = "sample_id")
plotReducedDim(sce_sub, dimred = "PCA", ncomponents = 2:3, colour_by = "sample_id")

# plot Vireo unassigned/doublets in PCA space: these have (almost all) been removed by filtering
plotReducedDim(sce_sub, dimred = "PCA", colour_by = "unassigned")
plotReducedDim(sce_sub, dimred = "PCA", ncomponents = 2:3, colour_by = "unassigned")

plotReducedDim(sce_sub, dimred = "PCA", colour_by = "doublet")
plotReducedDim(sce_sub, dimred = "PCA", ncomponents = 2:3, colour_by = "doublet")


# UMAP
plotReducedDim(sce_sub, dimred = "UMAP", colour_by = "sample_id")
plotReducedDim(sce_sub, dimred = "UMAP", colour_by = "unassigned")
plotReducedDim(sce_sub, dimred = "UMAP", colour_by = "doublet")
plotReducedDim(sce_sub, dimred = "UMAP", colour_by = "n_genes")
plotReducedDim(sce_sub, dimred = "UMAP", colour_by = "prop_mito")


# how many unassigned/doublets remain after filtering
sum((!colData(sce)$discard) & colData(sce)$unassigned)
sum((!colData(sce)$discard) & colData(sce)$doublet)



