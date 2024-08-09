library(tidyverse)
rm(list=ls())

# Set wd to parent experiment folder
setwd("~/Dropbox/projects/fmt/experiments/2021-03-18_JxE_FMT03/")

# Load Data
meta <- read_csv("metadata/samples.csv")
tax <- read_csv("data/tax.csv")
zotu <- read_csv("data/zotu.csv") 

# Make tidy
zotu <- zotu %>%
  pivot_longer(!seqid, names_to = "OTU", values_to = "read_count")

# Need to figure out spike-in strain, Sporosarcina, Otu id's
sporo.OTUs <- tax %>%
  filter(genus == "Sporosarcina") %>%
  select(OTU) %>%
  as_vector()

# Compute sporo values to normalize by for abs abd
sporo.norm <- zotu %>%
  filter(OTU %in% sporo.OTUs) %>%
  group_by(seqid) %>%
  summarise(sporo_read_count = sum(read_count))

# With sporo norm created, filter out sporo reads
# before calc'ing rel_abd
zotu <- zotu %>%
  filter(!OTU %in% sporo.OTUs)

# Add sporo norm values
zotu <- zotu %>%
  left_join(sporo.norm)

# Add pellet weights
zotu <- zotu %>%
  left_join(select(meta, seqid, pellet_weight))

# compute absolute abundance
# AU / g
# compute abs abd
zotu <- zotu %>%
  mutate(abs_abd = read_count / sporo_read_count / pellet_weight )

# Add relative abundance
zotu <- zotu %>%
  group_by(seqid) %>%
  mutate(rel_abd = read_count / sum(read_count)) %>%
  ungroup()

# sanity check
zotu %>%
  group_by(seqid) %>%
  summarise(total = sum(rel_abd)) %>%
  select(total) %>%
  unique()

# Pull out the taxonomy info I'm interested in
tax.to.add <- tax %>%
  select(OTU, domain, phylum, class, order, family, genus)  

zotu <- zotu %>%
  left_join(tax.to.add)

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


zotu %>%
  write_csv("data/tidy_zotu.csv")

