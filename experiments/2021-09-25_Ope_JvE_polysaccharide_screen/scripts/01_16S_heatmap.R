library(tidyverse)

# Clean up working space
rm(list=ls())

# Set working dir
setwd("~/Dropbox/projects/fmt/experiments/2021-09-25_Ope_JvE_polysaccharide_screen/")

# Load long form otu table
lotu <- read_csv("data/lotu.csv")

# Prep Work

# Pull out OTUs that are atleast 1% relative abundance
# in atleast one condition
OTU.filter <- lotu %>%
  group_by(mouse_type, polysacc, OTU) %>%
  summarise(mean_rel_abd = mean(rel_abd), .groups="drop") %>%
  filter(mean_rel_abd >= .01) %>%
  select(OTU) %>%
  unique() %>%
  as_vector()

# Filter to abundant OTUs
lotu <- lotu %>%
  filter(OTU %in% OTU.filter)

# Arrange phylogeny as a factor based on family
# for fixed order plotting
family.levels <-lotu %>%
  arrange(domain, phylum, class, order, family, genus) %>%
  select(family) %>%
  unique() %>%
  as_vector()

# convert family to factor
lotu$family <- factor(lotu$family, levels = family.levels)

# Add log10 data
lotu <- lotu %>%
  mutate(log10_rel_abd = if_else(rel_abd == 0, -6, log10(rel_abd)))

# Add mouse labels for plotting
lotu <- lotu %>%
  mutate(mouse = paste0(mouse_type, replicate))

# convert mouse labels into factor for fixed order plotting
lotu$mouse <- factor(lotu$mouse, levels = c("JC1", "JC2", "EC1", "EC2", "JE1", "JE2"))

# Create blue -> grey -> red pallette
pal.colors <- c("#647CBD", "#DBDDE0", "#E0624E")

# Plot
heatmap <- lotu %>%
  ggplot() +
  geom_tile(aes(x=mouse, y=OTU, fill=log10_rel_abd), color = "black") +
  facet_grid(family~polysacc, space="free", scale="free") +
  theme_bw(base_size = 8) + 
  theme(strip.text.y.right = element_text(angle = 0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradientn(colours = pal.colors) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

heatmap

ggsave("figures/heatmap.png", heatmap)
