library(tidyverse)
library(wesanderson)
rm(list=ls())

setwd("~/Dropbox/projects/fmt/experiments/2018-11-13_JxE_FMT02/")

meta <- read_csv("metadata/samples.csv")
tax <- read_csv("data/tax.csv")
zotu <- read_csv("data/zotu.csv")

# Make tidy
zotu <- zotu %>%
  pivot_longer(!seqid, names_to = "OTU", values_to = "read_count")

# Add relative abundance
zotu <- zotu %>%
  group_by(seqid) %>%
  mutate(rel_abd = read_count / sum(read_count)) %>%
  ungroup()

# Pull out the taxonomy info I'm interested in
tax.to.add <- tax %>%
  select(OTU, domain, phylum, class, order, family, genus)

zotu <- zotu %>%
  left_join(tax.to.add)

# sanity check
zotu %>%
  group_by(seqid) %>%
  summarise(total = sum(rel_abd)) %>%
  select(total) %>%
  unique()

# pull out useful metadata to add to zotu
meta.to.add <- meta %>%
  select(seqid, day, tx_group, mouse_type, mouse_num)

# add in metadata
zotu <- zotu %>%
  left_join(meta.to.add)

# Convert day to a factor, w/ levels to force order
day.levels <- zotu %>%
  select(day) %>%
  arrange(day) %>%
  unique() %>%
  as_vector()

zotu$day <- factor(zotu$day, levels = day.levels)

# Force phylo order
genus.levels <- zotu %>%
  arrange(domain, phylum, class, order, family, genus) %>%
  select(genus) %>%
  unique() %>%
  as_vector()

zotu$genus <- factor(zotu$genus, levels = genus.levels)

# Create list of OTUs filtered to certain average abundance at day 0
# done independently for both Jax & Envigo, then merged and unique'd
filter.level <- .01
filtered.OTUs <- zotu %>%
  filter(day == 0) %>%
  group_by(mouse_type, OTU) %>%
  summarise(mean_rel_abd = mean(rel_abd), .groups = "drop") %>%
  filter(mean_rel_abd > filter.level) %>%
  select(OTU) %>%
  unique() %>%
  as_vector()

### Visualization
# Prep color palletes
cpal <- wes_palette("Zissou1", 100, type = "continuous")

# Relative abundance heatmap : Day vs. OTU, color by log10 rel_abd
rel.abd.heatmap <- zotu %>%
  filter(OTU %in% filtered.OTUs) %>%
  filter(mouse_type == "Jackson") %>%
  group_by(tx_group, day, genus, OTU) %>%
  summarise(mean_rel_abd = mean(rel_abd), .groups = "keep") %>%
  ungroup() %>%
  mutate(log10_mean_rel_abd = if_else(mean_rel_abd == 0, -6, log10(mean_rel_abd))) %>%
  ggplot() + 
  geom_tile(aes(x=day, y=OTU, fill=log10_mean_rel_abd)) +
  facet_grid(genus~tx_group, space="free", scales="free") +
  scale_fill_gradientn(colours = cpal) +
  theme(strip.text.y.right = element_text(angle = 0))
rel.abd.heatmap
ggsave("figures/rel_abd_heatmap.svg", rel.abd.heatmap)
