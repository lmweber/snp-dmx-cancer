########################
# Estimated cost savings
########################

# Plot of estimated cost savings from demultiplexing

# Calculated using Satija Lab "Cost Per Cell" online calculator: 
# https://satijalab.org/costpercell


library(tidyverse)
library(magrittr)
library(ggplot2)


# ----------------------
# estimated cost savings
# ----------------------

# set up data frame using calculations from Satija Lab "Cost Per Cell" website

df_savings <- rbind(
  tibble(n_samples = 4, multiplexing = 3942, no_multiplexing = 9673), 
  tibble(n_samples = 6, multiplexing = 5459, no_multiplexing = 14523), 
  tibble(n_samples = 8, multiplexing = 8063, no_multiplexing = 19373)
)

df_plot <- gather(df_savings, "multiplexing", "cost", "multiplexing", "no_multiplexing")

df_plot$multiplexing %<>% factor(levels = c("no_multiplexing", "multiplexing"), 
                                 labels = c("no multiplexing", "multiplexing"))


# -----------
# create plot
# -----------

colors <- c("dodgerblue3", "firebrick1")

ggplot(df_plot, aes(x = n_samples, y = cost, 
                    group = multiplexing, shape = multiplexing, color = multiplexing)) + 
  geom_line() + 
  geom_point(size = 3) + 
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = c(4, 6, 8)) + 
  ylim(c(0, max(df_plot$cost))) + 
  xlab("number of samples") + 
  ylab("cost ($)") + 
  ggtitle("Estimated cost for experiment", 
          subtitle = "With and without multiplexing prior to library preparation") + 
  theme_bw()

ggsave("../../plots/estimated_cost_savings.pdf", width = 5, height = 3.75)
ggsave("../../plots/estimated_cost_savings.png", width = 5, height = 3.75)

