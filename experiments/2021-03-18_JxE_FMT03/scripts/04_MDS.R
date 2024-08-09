library(tidyverse)
rm(list=ls())

# Set wd to parent experiment folder
setwd("~/Dropbox/projects/fmt/experiments/2021-03-18_JxE_FMT03/")

# Load Data
zotu <- read_csv("data/tidy_zotu.csv")

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
  filter(tx_group %in% c("jax_control", "envigo_control")) %>%
  group_by(mouse_type, OTU) %>%
  summarise(mean_rel_abd = mean(rel_abd), .groups = "drop") %>%
  filter(mean_rel_abd > filter.level) %>%
  select(OTU) %>%
  unique() %>%
  as_vector()

filtered.zotu <- zotu %>%
  filter(OTU %in% filtered.OTUs)

filtered.zotu %>%
  select(rel_abd) %>%
  filter(rel_abd > 0) %>%
  min()

# Prepare palette colors
pal.colors <- c("#647CBD", "#DBDDE0", "#E0624E")

filtered.zotu %>%
  filter(tx_group == "jax_envigo_mix") %>%
  ggplot() +
  geom_tile(aes(x=seqid, y=OTU, fill=log10(rel_abd))) +
  scale_fill_gradientn(colours = pal.colors,
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  facet_grid(.~day, scales = "free", space = "free")
  

# Relative abundance heatmap : Day vs. OTU, color by log10 rel_abd
rel.abd.heatmap <- filtered.zotu %>%
  filter(!(mouse_type == "envigo" & tx_group == "jax_envigo_mix")) %>%
  group_by(tx_group, day, genus, OTU) %>%
  summarise(mean_rel_abd = mean(rel_abd), .groups = "keep") %>%
  ungroup() %>%
  mutate(log10_mean_rel_abd = if_else(mean_rel_abd == 0, -5, log10(mean_rel_abd))) %>%
  ggplot() + 
  geom_tile(aes(x=day, y=OTU, fill=log10_mean_rel_abd), color = "black") +
  facet_grid(genus~tx_group, space="free", scales="free") +
  theme_bw() +
  theme(strip.text.y.right = element_text(angle = 0)) +
  scale_fill_gradientn(colours = pal.colors,
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
rel.abd.heatmap
ggsave("figures/rel_abd_heatmap.svg", rel.abd.heatmap)
ggsave("figures/rel_abd_heatmap.png", rel.abd.heatmap)

filtered.zotu %>%
  filter(!(mouse_type == "jax" & tx_group == "jax_envigo_mix")) %>%
  filter(tx_group == "jax_envigo_mix") %>%
  select(mouse_type) %>%
  unique()

