library(tidyverse)

setwd("~/Dropbox/projects/fmt/experiments/2021-09-25_Ope_JvE_polysaccharide_screen/")

# Load Data
meta <- read_csv("metadata/metadata.csv")
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

# Add that taxonomy info
zotu <- zotu %>%
  left_join(tax.to.add)

# pull out useful metadata to add to zotu
meta.to.add <- meta %>%
  select(seqid, mouse_type, polysacc, replicate)

# add in metadata
zotu <- zotu %>%
  left_join(meta.to.add)

zotu %>%
  write_csv("data/lotu.csv")
