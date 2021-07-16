###############
# Runtime plots
###############

# Runtime plot for cellSNP step in benchmarking scenarios

# HGSOC dataset, 20pc doublets simulation


# module load conda_R/4.0
# Rscript runtimes.R


library(tidyverse)
library(magrittr)
library(ggplot2)


# --------------------------------------------------------
# load runtimes for cellSNP step in benchmarking scenarios
# --------------------------------------------------------

# steps with one value for all samples combined

fn <- "../../benchmarking/scenarios/HGSOC/20pc/1000GenomesFilt_cellSNPVireo/runtimes/runtime_1000GenomesFilt_cellSNPVireo_cellSNP_HGSOC_20pc.txt"
runtime_1000GenomesFilt_cellSNP <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))

fn <- "../../benchmarking/scenarios/HGSOC/20pc/1000GenomesUnfilt_cellSNPVireo/runtimes/runtime_1000GenomesUnfilt_cellSNPVireo_cellSNP_HGSOC_20pc.txt"
runtime_1000GenomesUnfilt_cellSNP <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))

fn <- "../../benchmarking/scenarios/HGSOC/20pc/bulkBcftools_cellSNPVireo/runtimes/runtime_bulkBcftools_cellSNPVireo_cellSNP_HGSOC_20pc.txt"
runtime_bulkBcftools_cellSNP <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))

fn <- "../../benchmarking/scenarios/HGSOC/20pc/bulkCellSNP_cellSNPVireo/runtimes/runtime_bulkCellSNP_cellSNPVireo_cellSNP_HGSOC_20pc.txt"
runtime_bulkCellSNP_cellSNP <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))

fn <- "../../benchmarking/scenarios/HGSOC/20pc/singlecellCellSNP_cellSNPVireo/runtimes/runtime_singlecellCellSNP_cellSNPVireo_cellSNP_HGSOC_20pc.txt"
runtime_singlecellCellSNP_cellSNP <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))


# -------------
# plot runtimes
# -------------

# set up plotting data frame

df_single <- as.data.frame(rbind(
  X1000GenomesFilt_cellSNP = runtime_1000GenomesFilt_cellSNP, 
  X1000GenomesUnfilt_cellSNP = runtime_1000GenomesUnfilt_cellSNP, 
  bulkBcftools_cellSNP = runtime_bulkBcftools_cellSNP, 
  bulkCellSNP_cellSNP = runtime_bulkCellSNP_cellSNP, 
  singlecellCellSNP_cellSNP = runtime_singlecellCellSNP_cellSNP
))
# fix names
rownames(df_single)[1:2] <- gsub("^X", "", rownames(df_single)[1:2])

colnames(df_single) <- "all"
df_single$all %<>% as.numeric
df_single$method <- rownames(df_single)
df_single$parallel <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
df_single <- gather(df_single, "sample_id", "runtime", "all")

df_combined <- df_single

method_names <- c("1000GenomesFilt_cellSNP", "1000GenomesUnfilt_cellSNP", "bulkBcftools_cellSNP", 
                  "bulkCellSNP_cellSNP", "singlecellCellSNP_cellSNP")

df_combined$method <- factor(df_combined$method, levels = method_names)
df_combined$sample_id <- factor(df_combined$sample_id, levels = c("all"))


# convert units to hours

df_plot <- df_combined
df_plot$runtime <- df_plot$runtime / 3600


# generate plot

colors <- c("firebrick")

ggplot(df_plot, aes(x = runtime, y = method, color = parallel)) + 
  geom_point(shape = 4, size = 1.5, stroke = 1.5) + 
  scale_color_manual(values = colors) + 
  xlim(c(0, max(df_plot$runtime))) + 
  xlab("runtime (hours)") + 
  scale_y_discrete(limits = rev(levels(df_plot$method))) + 
  ggtitle("cellSNP scenarios") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        legend.position = "none")

ggsave("../../plots/main/runtimes_cellSNP_scenarios_HGSOC_20pc.pdf", width = 4.5, height = 2.4)
ggsave("../../plots/main/runtimes_cellSNP_scenarios_HGSOC_20pc.png", width = 4.5, height = 2.4)

