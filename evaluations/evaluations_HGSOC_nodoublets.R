##################
# Evaluation plots
##################

# Script to generate plots for performance evaluations of benchmarking scenarios

# HGSOC dataset, no doublets


# module load conda_R/4.0
# Rscript evaluations.R


library(tidyverse)
library(magrittr)
library(ggplot2)


# -----------------
# load ground truth
# -----------------

# load ground truth stored in sample ID suffixes in cell barcodes

fn_barcodes_merged <- "../../benchmarking/outputs/HGSOC/barcodes_merged/barcodes_merged.tsv"
barcodes_merged <- read_table(fn_barcodes_merged, col_names = "barcode_id")

df_truth <- barcodes_merged %>% 
  mutate(sample_id = gsub("^[A-Za-z]+-", "", barcode_id))

head(df_truth)


# ----------------------
# benchmarking scenarios
# ----------------------

# Vireo scenarios: outputs saved in file "donor_ids.tsv"
# demuxlet scenarios: outputs saved in file

scenario_names <- c(
  "1000GenomesFilt_cellSNPVireo", 
  "1000GenomesUnfilt_cellSNPVireo", 
  "bulkBcftools_cellSNPVireo", 
  "bulkBcftools_demuxlet", 
  "bulkCellSNP_cellSNPVireo", 
  "singlecellCellSNP_cellSNPVireo"
)

summary_tables <- vector("list", length(scenario_names))
names(summary_tables) <- scenario_names

precision <- vector("list", length(scenario_names))
names(precision) <- scenario_names

recall <- vector("list", length(scenario_names))
names(recall) <- scenario_names


for (i in 1:length(scenario_names)) {
  
  # -----------
  # load output
  # -----------
  
  # Vireo scenarios
  if (i %in% c(1:3, 5:6)) {
    fn <- paste0("../../benchmarking/scenarios/HGSOC/nodoublets/", scenario_names[i], "/vireo/donor_ids.tsv")
    out <- read_tsv(fn)
    out_sub <- out[, c("cell", "donor_id")]
    
    stopifnot(nrow(df_truth) == nrow(out_sub))
    
    # match rows
    df_truth_tmp <- right_join(df_truth, out_sub, by = c("barcode_id" = "cell"))
    
    df_truth_tmp$truth <- factor(df_truth_tmp$sample_id)
    df_truth_tmp$predicted <- factor(df_truth_tmp$donor_id)
  }
  
  # demuxlet scenarios
  if (i == 4) {
    fn <- paste0("../../benchmarking/scenarios/HGSOC/nodoublets/", scenario_names[i], "/demuxlet.best")
    out <- read_tsv(fn)
    out_sub <- out[, c("BARCODE", "BEST")]
    
    stopifnot(nrow(df_truth) == nrow(out_sub))
    
    # match rows
    df_truth_tmp <- right_join(df_truth, out_sub, by = c("barcode_id" = "BARCODE"))
    
    df_truth_tmp$truth <- factor(df_truth_tmp$sample_id)
    df_truth_tmp$predicted <- factor(df_truth_tmp$BEST)
  }
  
  
  # --------------------------------------------
  # calculate summary table and match sample IDs
  # --------------------------------------------
  
  # summary table
  tbl_summary <- table(truth = df_truth_tmp$truth, predicted = df_truth_tmp$predicted)
  tbl_summary
  
  # manually match sample IDs (since predicted sample IDs are in arbitrary order)
  
  # Vireo scenarios
  if (i == 1) {
    levels(df_truth_tmp$predicted)[1:3] <- c("X3", "X4", "X2")
  } else if (i == 2) {
    levels(df_truth_tmp$predicted)[1:3] <- c("X2", "X4", "X3")
  } else if (i == 3) {
    levels(df_truth_tmp$predicted)[1:3] <- c("X2", "X4", "X3")
  } else if (i == 5) {
    levels(df_truth_tmp$predicted)[1:3] <- c("X2", "X3", "X4")
  } else if (i == 6) {
    levels(df_truth_tmp$predicted)[1:3] <- c("X2", "X3", "X4")
  }
  
  # demuxlet scenarios
  if (i == 4) {
    levels(df_truth_tmp$predicted)[9:11] <- c("X3", "X2", "X4")
  }
  
  # updated summary table
  tbl_summary <- table(truth = df_truth_tmp$truth, predicted = df_truth_tmp$predicted)
  tbl_summary
  
  summary_tables[[i]] <- tbl_summary
  
  
  # ------------------------------
  # calculate precision and recall
  # ------------------------------
  
  # calculate recall
  re <- sapply(rownames(tbl_summary), function(s) {
    tbl_summary[s, s] / sum(tbl_summary[s, ])
  })
  
  # calculate precision
  pr <- sapply(rownames(tbl_summary), function(s) {
    tbl_summary[s, s] / sum(tbl_summary[, s])
  })
  
  recall[[i]] <- re
  precision[[i]] <- pr
  
}


# ------------------------------
# format data frame for plotting
# ------------------------------

df_recall <- as.data.frame(t(data.frame(recall, check.names = FALSE)))
df_recall$scenario <- rownames(df_recall)
df_recall$metric <- "recall"

df_precision <- as.data.frame(t(data.frame(precision, check.names = FALSE)))
df_precision$scenario <- rownames(df_precision)
df_precision$metric <- "precision"

df_plot <- rbind(df_recall, df_precision)

df_plot <- gather(df_plot, "sample_id", "value", "X2", "X3", "X4")

df_plot$scenario <- factor(df_plot$scenario)
df_plot$metric <- as.factor(df_plot$metric)
df_plot$sample_id <- as.factor(df_plot$sample_id)

df_plot <- spread(df_plot, "metric", "value")


# --------------
# generate plots
# --------------

# color palette (modified Okabe-Ito)
pal <- unname(palette.colors(palette = "Okabe-Ito"))
pal[1] <- "darkmagenta"

ggplot(df_plot, aes(x = recall, y = precision, color = scenario, shape = sample_id)) + 
  geom_point(size = 1.5, stroke = 1) + 
  scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1, 2, 0)) + 
  ggtitle("Precision-recall: HGSOC, no doublets") + 
  theme_bw()

ggsave("../../plots/precision_recall_HGSOC_nodoublets.pdf", width = 6.25, height = 3.5)
ggsave("../../plots/precision_recall_HGSOC_nodoublets.png", width = 6.25, height = 3.5)

