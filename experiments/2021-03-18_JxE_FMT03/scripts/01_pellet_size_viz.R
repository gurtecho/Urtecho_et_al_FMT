library(tidyverse)

# Set working dir to the parent experiment folder
setwd("~/Dropbox/projects/fmt/experiments/2021-03-18_JxE_replication/")

# Load the samples data sheet that contains pellet weights
weights <- read_csv("metadata/samples.csv")

# Visualize pellet size distributions
weights %>%
  filter(!is.na(pellet_weight)) %>%
  ggplot() +
  geom_boxplot(aes(x=tx_group, y=pellet_weight))
