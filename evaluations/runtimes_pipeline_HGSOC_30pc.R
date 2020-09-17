###############
# Runtime plots
###############

# Runtime plots for steps in main pipeline (note: not including bulk sample
# genotyping, which is shown in separate plot)

# HGSOC dataset, 30pc doublets simulation, best-performing scenario


# module load conda_R/4.0
# Rscript evaluations.R


library(tidyverse)
library(magrittr)
library(ggplot2)


# --------------------------------
# load runtimes for pipeline steps
# --------------------------------

# initial steps: steps with one value per sample

sample_names <- c("X2", "X3", "X4")
sample_names_long <- c("16030X2_HJVMLDMXX", "16030X3_HJTWLDMXX", "16030X4_HJTWLDMXX")

runtime_cellranger <- vector("list", 3)
names(runtime_cellranger) <- sample_names

runtime_parse_BAM_files <- vector("list", 3)
names(runtime_parse_BAM_files) <- sample_names

for (i in 1:length(sample_names)){
  fn <- paste0("../../benchmarking/runtimes/HGSOC/cellranger/runtime_cellranger_", sample_names_long[i], ".txt")
  runtime_cellranger[[i]] <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))
  
  fn <- paste0("../../benchmarking/runtimes/HGSOC/parse_BAM_files/runtime_parse_BAM_files_", sample_names_long[i], ".txt")
  runtime_parse_BAM_files[[i]] <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))
}

runtime_cellranger <- unlist(runtime_cellranger)
runtime_parse_BAM_files <- unlist(runtime_parse_BAM_files)


# initial steps: steps with one value for all samples combined

fn <- "../../benchmarking/runtimes/HGSOC/merge_and_index_BAM/runtime_merge_and_index_BAM.txt"
runtime_merge_and_index_BAM <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))

fn <- "../../benchmarking/runtimes/HGSOC/parse_and_merge_barcodes/runtime_parse_and_merge_barcodes.txt"
runtime_parse_and_merge_barcodes <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))


# doublets simulation

fn <- "../../benchmarking/scenarios/HGSOC/30pc/runtime_lookup_table_doublets_HGSOC_30pc.txt"
runtime_lookup_table_doublets <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))

fn <- "../../benchmarking/scenarios/HGSOC/30pc/runtime_parse_BAM_doublets_HGSOC_30pc.txt"
runtime_parse_BAM_doublets <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))


# demultiplexing: best-performing scenario (bulkBcftools_cellSNPVireo)

fn <- "../../benchmarking/scenarios/HGSOC/30pc/bulkBcftools_cellSNPVireo/runtimes/runtime_bulkBcftools_cellSNPVireo_cellSNP_HGSOC_30pc.txt"
runtime_cellSNP <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))

fn <- "../../benchmarking/scenarios/HGSOC/30pc/bulkBcftools_cellSNPVireo/runtimes/runtime_bulkBcftools_cellSNPVireo_vireo_HGSOC_30pc.txt"
runtime_vireo <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))


# --------------------------------
# plot runtimes for pipeline steps
# --------------------------------

# set up plotting data frame

df_separate <- as.data.frame(rbind(
  cellranger = runtime_cellranger, 
  parse_BAM_files = runtime_parse_BAM_files
))
for (i in 1:ncol(df_separate)) {
  df_separate[, i] %<>% as.numeric
}
df_separate$method <- rownames(df_separate)
df_separate <- gather(df_separate, "sample_id", "runtime", "X2", "X3", "X4")

df_single <- as.data.frame(rbind(
  merge_and_index_BAM = runtime_merge_and_index_BAM, 
  parse_and_merge_barcodes = runtime_parse_and_merge_barcodes, 
  lookup_table_doublets = runtime_lookup_table_doublets, 
  parse_BAM_doublets = runtime_parse_BAM_doublets, 
  cellSNP = runtime_cellSNP, 
  vireo = runtime_vireo
))
colnames(df_single) <- "all"
df_single$all %<>% as.numeric
df_single$method <- rownames(df_single)
df_single <- gather(df_single, "sample_id", "runtime", "all")

df_combined <- rbind(df_separate, df_single)

method_names <- c("cellranger", "parse_BAM_files", "merge_and_index_BAM", "parse_and_merge_barcodes", 
                  "lookup_table_doublets", "parse_BAM_doublets", "cellSNP", "vireo")

df_combined$method <- factor(df_combined$method, levels = method_names)
df_combined$sample_id <- factor(df_combined$sample_id, levels = c("X2", "X3", "X4", "all"))


# convert units to hours

df_plot <- df_combined
df_plot$runtime <- df_plot$runtime / 3600


# generate plot

ggplot(df_plot, aes(x = method, y = runtime, group = sample_id)) + 
  geom_point(color = "#D55E00", shape = 4, size = 1.5, stroke = 1.5) + 
  ylab("runtime (hours)") + 
  ggtitle("Runtimes: pipeline") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../../plots/runtimes_pipeline_HGSOC_30pc.pdf", width = 3.5, height = 4.5)
ggsave("../../plots/runtimes_pipeline_HGSOC_30pc.png", width = 3.5, height = 4.5)

