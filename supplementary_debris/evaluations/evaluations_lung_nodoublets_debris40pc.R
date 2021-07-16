##################
# Evaluation plots
##################

# Script to generate plots for performance evaluations of benchmarking scenarios

# lung dataset, no doublets / 40pc debris simulation


# module load conda_R/4.0
# Rscript evaluations.R


library(tidyverse)
library(magrittr)
library(ggplot2)


# -----------------
# load ground truth
# -----------------

# load ground truth stored in sample ID suffixes in cell barcodes

# barcodes
fn_barcodes_merged <- "../../../supplementary_debris/scenarios/lung/nodoublets/barcodes_merged_lung_nodoublets_debris40pc.tsv"
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
  "1000GenomesFilt_cellSNPVireo"
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
  if (i == 1) {
    fn <- paste0("../../../supplementary_debris/scenarios/lung/nodoublets/debris40pc/", scenario_names[i], "/vireo/donor_ids.tsv")
    out <- read_tsv(fn)
    out_sub <- out[, c("cell", "donor_id")]
    
    stopifnot(nrow(df_truth) == nrow(out_sub))
    
    # match rows
    df_truth_tmp <- right_join(df_truth, out_sub, by = c("barcode_id" = "cell"))
    
    df_truth_tmp$truth <- factor(df_truth_tmp$sample_id)
    df_truth_tmp$predicted <- factor(df_truth_tmp$donor_id)
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
    levels(df_truth_tmp$predicted)[1:6] <- c("T28", "T31", "T20", "T08", "T25", "T09")
  }
  
  # updated summary table
  tbl_summary <- table(truth = df_truth_tmp$truth, predicted = df_truth_tmp$predicted)
  tbl_summary
  
  summary_tables[[i]] <- tbl_summary
  
  
  # ------------------------------
  # calculate precision and recall
  # ------------------------------
  
  # calculate recall
  re <- sapply(rownames(tbl_summary)[1:6], function(s) {
    tbl_summary[s, s] / sum(tbl_summary[s, ])
  })
  
  # calculate precision
  pr <- sapply(rownames(tbl_summary)[1:6], function(s) {
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

df_plot <- gather(df_plot, "sample_id", "value", "T08", "T09", "T20", "T25", "T28", "T31")

df_plot$scenario <- factor(df_plot$scenario)
df_plot$metric <- as.factor(df_plot$metric)
df_plot$sample_id <- as.factor(df_plot$sample_id)

# summary values
df_plot %>% 
  group_by(scenario, metric) %>% 
  summarize(mean = mean(value))

# data frame for plotting
df_plot <- spread(df_plot, "metric", "value")


# --------------
# generate plots
# --------------

# color palette (modified Okabe-Ito)
pal <- unname(palette.colors(palette = "Okabe-Ito"))
pal[1] <- "darkmagenta"
pal <- pal[1]

ggplot(df_plot, aes(x = recall, y = precision, color = scenario, shape = sample_id)) + 
  geom_point(size = 1.5, stroke = 1) + 
  scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(1, 2, 0, 3, 4, 5)) + 
  xlim(0, 1) + 
  ylim(0, 1) + 
  ggtitle("Lung, no doublets, 40% debris") + 
  theme_bw()

ggsave("../../../plots/supp_debris/precision_recall_lung_nodoublets_debris40pc.pdf", width = 6.15, height = 3.5)
ggsave("../../../plots/supp_debris/precision_recall_lung_nodoublets_debris40pc.png", width = 6.15, height = 3.5)

