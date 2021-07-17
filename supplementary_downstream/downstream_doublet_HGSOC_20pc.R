#####################################################################
# Script to run downstream doublet detection and evaluate performance
# Lukas Weber, July 2021
#####################################################################

# Script to run downstream doublet detection tool and evaluate performance for 
# one simulation scenario


library(DropletUtils)
library(readr)
library(tidyverse)
library(scater)
library(scran)
library(BiocSingular)


# ---------
# load data
# ---------

# load data and create SingleCellExperiment

samples <- c(
  "../../benchmarking/outputs/HGSOC/16030X2_HJVMLDMXX/outs/filtered_feature_bc_matrix", 
  "../../benchmarking/outputs/HGSOC/16030X3_HJTWLDMXX/outs/filtered_feature_bc_matrix", 
  "../../benchmarking/outputs/HGSOC/16030X4_HJTWLDMXX/outs/filtered_feature_bc_matrix"
)
sample_names <- c("16030X2", "16030X3", "16030X4")

sce <- read10xCounts(samples, sample.names = sample_names)

# update cell barcodes to include sample IDs to match format in merged BAM files

colData(sce)$sample_id <- gsub("([0-9]+)(.*)", "\\2", colData(sce)$Sample)
colData(sce)$barcode_id <- paste0(gsub("([A-Z]+-)1", "\\1", colData(sce)$Barcode), colData(sce)$sample_id)


# -----------------
# identify doublets
# -----------------

# identify true doublets

fn_lookup_table <- "../../benchmarking/scenarios/HGSOC/20pc/lookup_table_doublets_HGSOC_20pc.tsv"
lookup_table <- read_tsv(fn_lookup_table)

# remove barcodes that were merged into other barcodes when generating doublets
sce <- sce[, !(colData(sce)$barcode_id %in% lookup_table$original)]
dim(sce)

# identify doublets
colData(sce)$true_doublet <- colData(sce)$barcode_id %in% lookup_table$replacement


# identify doublets detected by Vireo

# for best-performing scenario: bulkBcftools_cellSNPVireo

# Vireo outputs are saved in file "donor_ids.tsv"
fn <- "../../benchmarking/scenarios/HGSOC/20pc/bulkBcftools_cellSNPVireo/vireo/donor_ids.tsv"
out_vireo <- read_tsv(fn)

barcodes_dbl_vireo <- out_vireo$cell[out_vireo$donor_id == "doublet"]

colData(sce)$vireo_doublet <- colData(sce)$barcode_id %in% barcodes_dbl_vireo


# subset to keep cells that were not identified as doublets by Vireo, since these 
# are the cells where we want to apply downstream doublet detection methods

table(true = colData(sce)$true_doublet, vireo = colData(sce)$vireo_doublet)
table(colData(sce)$vireo_doublet)
dim(sce)

sce <- sce[, !colData(sce)$vireo_doublet]
dim(sce)


# ----------------------------
# preprocessing and clustering
# ----------------------------

# run preprocessing steps and clustering
# from https://bioconductor.org/books/release/OSCA/doublet-detection.html

# quality control

is_mito <- grepl("^MT-", rowData(sce)$Symbol)
table(is_mito)

# note: skip quality control filtering since all cells were used for generating doublets

# stats <- perCellQCMetrics(sce, subsets = list(mito = which(rowData(sce)$is_mito)))
# qc <- quickPerCellQC(stats, percent_subsets = "subsets_mito_percent")
# sce <- sce[, !qc$discard]

# normalization

set.seed(123)
qclus <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = qclus)
sce <- logNormCounts(sce)

# feature selection

sce <- sce[!is_mito, ]
dim(sce)

set.seed(123)
dec <- modelGeneVar(sce)
top <- getTopHVGs(dec, prop = 0.1)

# dimensionality reduction

set.seed(123)
sce <- runPCA(sce, subset_row = top)
#sce <- runTSNE(sce, dimred = "PCA")

# clustering

set.seed(123)
snn_gr <- buildSNNGraph(sce, use.dimred = "PCA", k = 10)
colLabels(sce) <- factor(igraph::cluster_walktrap(snn_gr)$membership)
table(colLabels(sce))


# --------------------------------
# run downstream doublet detection
# --------------------------------

# see https://bioconductor.org/books/release/OSCA/doublet-detection.html

library(scDblFinder)

# run doublet detection method

dbl_out <- findDoubletClusters(sce)
dbl_out

# select clusters likely to represent doublets using elbow method by plotting num.de

plot(dbl_out$num.de)

threshold <- 70


# formatted plot
df <- data.frame(
  cluster = rownames(dbl_out), 
  num.de = dbl_out$num.de
)
ggplot(df, aes(x = cluster, y = num.de)) + 
  geom_point() + 
  scale_x_discrete(name = "cluster ID", limits = df$cluster) + 
  geom_hline(yintercept = threshold, color = "red") + 
  labs(y = "number of DE genes") + 
  ggtitle("HGSOC, 20% doublets") + 
  theme_bw()

ggsave("../../plots/supp_downstream/downstream_doublets_elbow_HGSOC_20pc.pdf", width = 7, height = 4.5)
ggsave("../../plots/supp_downstream/downstream_doublets_elbow_HGSOC_20pc.png", width = 7, height = 4.5)


# choose elbow
chosen_dbl <- rownames(dbl_out[dbl_out$num.de < threshold, ])
chosen_dbl

# alternatively: using outlier method
# chosen_dbl <- rownames(dbl_out)[isOutlier(dbl_out$num.de, type = "lower", log = TRUE)]
# chosen_dbl

colData(sce)$chosen_doublet <- colData(sce)$label %in% chosen_dbl

table(true = colData(sce)$true_doublet, chosen = colData(sce)$chosen_doublet)

# output:
#        chosen
# true   FALSE TRUE
# FALSE  6405  3720
# TRUE   1003   600

