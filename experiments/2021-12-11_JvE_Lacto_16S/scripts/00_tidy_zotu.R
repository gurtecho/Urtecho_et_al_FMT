library(tidyverse)
library(RColorBrewer)

# Clean up env
rm(list =ls())

# set wd to the experiment folder
setwd("~/Dropbox/projects/fmt/experiments/2021-12-11_JvE_Lacto_16S/")

# Load data sets
zotu <- read_csv("data/zotu.csv")
tax <- read_csv("data/tax.csv")
meta <- read_csv("metadata/metadata.csv")

# Covert zotu into tidy format
tidy.zotu <- zotu %>%
  pivot_longer(-seqid, names_to = "ZOTU", values_to = "count")

# Pull out list of spike in ZOTUs
spike.in <- tax %>%
  filter(genus == "Sporosarcina") %>%
  select(OTU) %>%
  as_vector()

# Pull out spike-in strains and generate table with counts
spike.count <- tidy.zotu %>%
  filter(ZOTU %in% spike.in) %>%
  group_by(seqid) %>%
  summarise(spike_count = sum(count))

# Filter spike-in counts from the ZOTU table
tidy.zotu <- tidy.zotu %>%
  filter(!ZOTU %in% spike.in)

# Compute relative abundance
tidy.zotu <- tidy.zotu %>%
  group_by(seqid) %>%
  mutate(rel_ab = count / sum(count))

# Sanity check
tidy.zotu %>%
  group_by(seqid) %>%
  summarise(total_ab = sum(rel_ab))

# Add in sporo count
tidy.zotu <- tidy.zotu %>%
  left_join(spike.count)

# filter the tax table and limit to desired tax levels
filt.tax <- tax %>% 
  select(OTU, order, family, genus) %>%
  rename(ZOTU = OTU)


# filter the metadata to desired columns
filt.meta <- meta %>%
  select(-c("S3_R1_path" ,"S3_R2_path"))

tidy.zotu <- tidy.zotu %>%
  left_join(filt.meta) %>%
  left_join(filt.tax)

filt.tidy.zotu <- tidy.zotu %>%
  filter(rel_ab >= .01)

order.barplot <- filt.tidy.zotu %>%
  ggplot() +
  geom_bar(aes(x = seqid, y = rel_ab, fill = order), stat = "identity") +
  facet_grid(.~mouse_type, space = "free", scale = "free")
order.barplot
ggsave("./figures/order_barplot.png", order.barplot)

lacto.zotu.barplot <- filt.tidy.zotu %>%
  filter(family == "Lactobacillaceae") %>%
  ggplot() +
  geom_bar(aes(x = seqid, y = rel_ab, fill = ZOTU), stat = "identity") +
  facet_grid(.~mouse_type, space = "free", scale = "free") +
  scale_fill_brewer(palette = "Accent") +
  ggtitle("Relative Abundance of Lactobacillaceae ZOTUs")
lacto.zotu.barplot
ggsave("./figures/lacto_zotu_barplot.png", lacto.zotu.barplot)

filt.tidy.zotu %>%
  filter(family == "Lactobacillaceae") %>%
  write_csv("metadata/lacto_zotus.csv")


