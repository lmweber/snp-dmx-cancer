###############
# Runtime plots
###############

# Runtime plot for bulk sample alignment steps

# HGSOC dataset


# module load conda_R/4.1
# Rscript runtimes.R


library(tidyverse)
library(magrittr)
library(ggplot2)


# ---------------------------------------------
# load runtimes for bulk sample alignment steps
# ---------------------------------------------

# steps with one value per sample

sample_names_bulk <- c("17667X1", "17667X2", "17667X3")

runtime_align_index_bulk <- vector("list", 3)
names(runtime_align_index_bulk) <- sample_names_bulk

for (i in 1:length(sample_names_bulk)){
  fn <- paste0("../../genotype/runtimes/align_index_bulk_STAR/runtime_align_index_bulk_", sample_names_bulk[i], ".txt")
  runtime_align_index_bulk[[i]] <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))
}

runtime_align_index_bulk <- unlist(runtime_align_index_bulk)


# steps with one value for all samples combined

fn <- "../../genotype/runtimes/align_index_bulk_STAR/runtime_create_STAR_index.txt"
runtime_create_STAR_index <- gsub("runtime: ", "", gsub(" seconds", "", readLines(fn)))


# -------------
# plot runtimes
# -------------

# set up plotting data frame

df_separate_bulk <- as.data.frame(rbind(
  align_index_bulk = runtime_align_index_bulk
))
for (i in 1:ncol(df_separate_bulk)) {
  df_separate_bulk[, i] %<>% as.numeric
}
df_separate_bulk$method <- rownames(df_separate_bulk)
df_separate_bulk$parallel <- c(TRUE)
df_separate_bulk <- gather(df_separate_bulk, "sample_id", "runtime", 
                           "17667X1", "17667X2", "17667X3")

df_single <- as.data.frame(rbind(
  create_STAR_index = runtime_create_STAR_index
))
colnames(df_single) <- "all"
df_single$all %<>% as.numeric
df_single$method <- rownames(df_single)
df_single$parallel <- c(TRUE)
df_single <- gather(df_single, "sample_id", "runtime", "all")

df_combined <- rbind(df_separate_bulk, df_single)

method_names <- c("create_STAR_index", "align_index_bulk")

df_combined$method <- factor(df_combined$method, levels = method_names)
df_combined$sample_id <- factor(df_combined$sample_id, 
                                levels = c("17667X1", "17667X2", "17667X3", "all"))


# convert units to hours

df_plot <- df_combined
df_plot$runtime <- df_plot$runtime / 3600


# generate plot

colors <- c("firebrick")

ggplot(df_plot, aes(x = runtime, y = method, group = sample_id, color = parallel)) + 
  geom_point(shape = 4, size = 1.5, stroke = 1.5) + 
  scale_color_manual(values = colors) + 
  xlim(c(0, max(df_plot$runtime))) + 
  xlab("runtime (hours)") + 
  scale_y_discrete(limits = rev(levels(df_plot$method))) + 
  ggtitle("Bulk samples alignment") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        legend.position = "none")

ggsave("../../plots/main/runtimes_bulk_samples_HGSOC.pdf", width = 4, height = 1.5)
ggsave("../../plots/main/runtimes_bulk_samples_HGSOC.png", width = 4, height = 1.5)

