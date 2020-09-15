###############
# Runtime plots
###############

# Runtime plots for genotyping tools

# HGSOC dataset


# module load conda_R/4.0
# Rscript evaluations.R


library(tidyverse)
library(magrittr)
library(ggplot2)


# ----------------------------------
# load runtimes for genotyping tools
# ----------------------------------

# steps with one value per sample

sample_names_bulk <- c("17667X1", "17667X2", "17667X3")
sample_names_singlecell <- c("16030X2", "16030X3", "16030X4")

runtime_align_index_bulk <- vector("list", 3)
names(runtime_align_index_bulk) <- sample_names_bulk

runtime_genotype_bulk_cellSNP <- vector("list", 3)
names(runtime_genotype_bulk_cellSNP) <- sample_names_bulk

runtime_genotype_singlecell_cellSNP <- vector("list", 3)
names(runtime_genotype_singlecell_cellSNP) <- sample_names_singlecell

for (i in 1:length(sample_names_bulk)){
  fn <- paste0("../../genotype/runtimes/align_index_bulk_STAR/runtime_align_index_bulk_", sample_names_bulk[i], ".txt")
  runtime_align_index_bulk[[i]] <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))
  
  fn <- paste0("../../genotype/runtimes/genotype_bulk_cellSNP/runtime_genotype_bulk_cellSNP_", sample_names_bulk[i], ".txt")
  runtime_genotype_bulk_cellSNP[[i]] <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))
  
  fn <- paste0("../../genotype/runtimes/genotype_singlecell_cellSNP/runtime_genotype_singlecell_cellSNP_", sample_names_singlecell[i], ".txt")
  runtime_genotype_singlecell_cellSNP[[i]] <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))
}

runtime_align_index_bulk <- unlist(runtime_align_index_bulk)
runtime_genotype_bulk_cellSNP <- unlist(runtime_genotype_bulk_cellSNP)
runtime_genotype_singlecell_cellSNP <- unlist(runtime_genotype_singlecell_cellSNP)

# steps with one value for all samples combined

fn <- "../../genotype/runtimes/align_index_bulk_STAR/runtime_create_STAR_index.txt"
runtime_create_STAR_index <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))

fn <- "../../genotype/runtimes/genotype_bulk_bcftools/runtime_genotype_bulk_HGSOC_bcftools.txt"
runtime_genotype_bulk_bcftools <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))


# ----------------------------------
# plot runtimes for genotyping tools
# ----------------------------------

# set up plotting data frame

df_separate_bulk <- as.data.frame(rbind(
  align_index_bulk = runtime_align_index_bulk, 
  genotype_bulk_cellSNP = runtime_genotype_bulk_cellSNP
))
for (i in 1:ncol(df_separate_bulk)) {
  df_separate_bulk[, i] %<>% as.numeric
}
df_separate_bulk$method <- rownames(df_separate_bulk)
df_separate_bulk <- gather(df_separate_bulk, "sample_id", "runtime", 
                           "17667X1", "17667X2", "17667X3")

df_separate_singlecell <- as.data.frame(rbind(
  genotype_singlecell_cellSNP = runtime_genotype_singlecell_cellSNP
))
for (i in 1:ncol(df_separate_singlecell)) {
  df_separate_singlecell[, i] %<>% as.numeric
}
df_separate_singlecell$method <- rownames(df_separate_singlecell)
df_separate_singlecell <- gather(df_separate_singlecell, "sample_id", "runtime", 
                                 "16030X2", "16030X3", "16030X4")

df_single <- as.data.frame(rbind(
  create_STAR_index = runtime_create_STAR_index, 
  genotype_bulk_bcftools = runtime_genotype_bulk_bcftools
))
colnames(df_single) <- "all"
df_single$all %<>% as.numeric
df_single$method <- rownames(df_single)
df_single <- gather(df_single, "sample_id", "runtime", "all")

df_combined <- rbind(df_separate_bulk, df_separate_singlecell, df_single)

method_names <- c("create_STAR_index", "align_index_bulk", "genotype_bulk_cellSNP", 
                  "genotype_bulk_bcftools", "genotype_singlecell_cellSNP")

df_combined$method <- factor(df_combined$method, levels = method_names)
df_combined$sample_id <- factor(df_combined$sample_id, 
                                levels = c("17667X1", "17667X2", "17667X3", 
                                           "16030X2", "16030X3", "16030X4", "all"))


# generate plot

ggplot(df_combined, aes(x = method, y = runtime, group = sample_id)) + 
  geom_point(color = "#D55E00", shape = 1, size = 1.5, stroke = 1.5) + 
  scale_y_log10() + 
  ylab("runtime (seconds)") + 
  ggtitle("Runtimes: genotyping tools") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../../plots/runtimes_genotype_HGSOC.pdf", width = 3, height = 4.5)
ggsave("../../plots/runtimes_genotype_HGSOC.png", width = 3, height = 4.5)

