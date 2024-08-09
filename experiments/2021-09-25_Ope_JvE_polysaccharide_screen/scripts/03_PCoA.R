library(tidyverse)
library(vegan)

rm(list=ls())

# Set wd to parent experiment folder
setwd("~/Dropbox/projects/fmt/experiments/2021-09-25_Ope_JvE_polysaccharide_screen/")

meta <- read_csv("metadata/metadata.csv")
lotu <- read_csv("data/lotu.csv")

community.matrix <- lotu %>%
  pivot_wider(id_cols = "seqid", names_from = "OTU", values_from = "rel_abd")

distance.matrix <- community.matrix %>%
  column_to_rownames(var = "seqid") %>%
  vegdist()

pcoa <- distance.matrix %>%
  cmdscale(k = 2, eig = TRUE)

pcoa.points <- pcoa$points %>%
  as_tibble(rownames = "seqid") %>%
  left_join(meta)

pcoa.points %>%
  ggplot() +
  geom_point(aes(x = V1, y = V2, color=mouse_type)) +
  theme_bw()

pcoa$eig[1] / sum(pcoa$eig)
pcoa$eig[2] / sum(pcoa$eig)
