library(tidyverse)
library("RColorBrewer")
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
filter.level <- .005
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

J <- filtered.zotu %>%
  filter(tx_group == "jax_control")
  

E <- filtered.zotu %>%
  filter(tx_group == "envigo_control")

JE <- filtered.zotu %>%
  filter(tx_group == "jax_envigo_mix") %>%
  filter(mouse_type == "jax")

exp.set <- J %>%
  bind_rows(E) %>%
  bind_rows(JE) %>%
  mutate(log10_rel_abd = if_else(rel_abd == 0, -5, log10(rel_abd)))
  
exp.matrix <- exp.set %>%
  pivot_wider(id_cols = "seqid", names_from = "OTU", values_from = "log10_rel_abd") %>%
  column_to_rownames(var = "seqid") %>%
  as.matrix()

# Prepare palette colors
pal.colors <- c("#647CBD", "#DBDDE0", "#E0624E")
pal.colors <- c("#647CBD", "#E0624E")
pal <- colorRampPalette(pal.colors)(256)

heatmap(t(exp.matrix), col = pal)

# Relative abundance heatmap : Day vs. OTU, color by log10 rel_abd
rel.abd.heatmap <- avg %>%
  ggplot() + 
  geom_tile(aes(x=day, y=OTU, fill=log10_mean_rel_abd), color = "black") +
  facet_grid(.~tx_group, space="free", scales="free") +
  theme_bw() +
  ylab("ZOTU") +
  xlab("Day post-Co-housing") +
  scale_fill_gradientn(colours = pal.colors,
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme(
    strip.text.y.right = element_text(angle = 0),
    strip.background=element_blank(), 
    panel.spacing = unit(0, "mm"), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
  )
rel.abd.heatmap
ggsave("figures/02_rel_abd_heatmap.svg", rel.abd.heatmap)
ggsave("figures/02_rel_abd_heatmap.png", rel.abd.heatmap)

filtered.zotu %>%
  filter(!(mouse_type == "jax" & tx_group == "jax_envigo_mix")) %>%
  filter(tx_group == "jax_envigo_mix") %>%
  select(mouse_type) %>%
  unique()

