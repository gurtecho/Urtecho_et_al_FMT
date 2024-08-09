library(tidyverse)

rm(list=ls())

# Working directory for the experiment
setwd("~/Dropbox/projects/fmt/experiments/2021-04-06_polysacc_growth_kinetics/")

#  Load the metadata
meta <- read_csv("metadata/metadata.csv")

# Load the tidy OD600 data processed w/ python script
data <- read_csv("data/data.csv")

# filter the data to those wells in the metadata
data <- data %>%
  filter(well %in% meta$well)

# Add the metadata to each entry
data <- data %>%
  left_join(meta, by="well")

# Compute the bg from the OD600 of the water wells
bg <- data %>%
  filter(polysacc == "water") %>%
  group_by(time_mins, sample_type, sample_dilution) %>%
  summarise(bg = mean(OD600))

# Add bg data to the actual data
# Remove water rows 
data <- data %>%
  left_join(bg, by=c("time_mins", "sample_type", "sample_dilution"))

# Compute OD600 normalized against water controls
data$norm_OD600 <- data$OD600 / data$bg

# Compute the OD600 corrected against water controls
data$corr_OD600 <- data$OD600 - data$bg

# Change the sample type to pretty
# version for plotting
# Make tx_group pretty
data <- data %>%
  mutate(sample_type = case_when(sample_type == "jax_control" ~ "Jax", 
                                sample_type == "envigo_control" ~ "Envigo"))


# convert sample type to factor which levels fixed in order I want to plot
data$sample_type = factor(data$sample_type, levels = c("Jax", "Envigo"))

# Make the polysacc pretty (capitalize)
data <- data %>%
  mutate(polysacc = if_else(polysacc == "mcf", "MCF", str_to_title(polysacc)))

# Order I want to plot polysacc in, to be used as levels
polysacc.order <-c("Water", "Arabinan", "Arabinogalactan",  "Arabinose", "Cellobiose", "Dextran",
                   "Inulin", "Pectic Galactan", "Pectin", "MCF")

# convert polysacc to factor, levels in order I want to plot
data$polysacc <- factor(data$polysacc, levels = polysacc.order)

# Convert Dilutions to pretty version
data <- data %>%
  mutate(sample_dilution = case_when(sample_dilution == 1.00 ~ "10^0", 
                                     sample_dilution == 0.10 ~ "10^-1",
                                     sample_dilution == 0.01 ~ "10^-2"))

# fix sample dilution order with a level'd factor
data$sample_dilution = factor(data$sample_dilution, levels = c("10^0", "10^-1", "10^-2"))

# Plot raw data
raw.viz <- data  %>%
  ggplot() +
  geom_line(aes(x=time_mins, y=OD600, color=polysacc)) +
  facet_grid(sample_type~sample_dilution) +
  xlab("Time (mins)") +
  labs(color="Polysaccharide") +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size=16),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size=16),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16),
    strip.text.x = element_text(size=16),
    strip.text.y = element_text(size=16))

raw.viz
ggsave("./figures/raw_data_lineplot.svg", raw.viz)
  

# Plot normalized data
norm.viz <- data  %>%
  filter(polysacc != "water") %>%
  ggplot() +
  geom_line(aes(x=time_mins, y=norm_OD600, color=polysacc)) +
  facet_grid(sample_type~sample_dilution) +
  xlab("Time (mins)") +
  ylab("Normalized OD600") +
  labs(color="Polysaccharide") +
  ylim(0, 4) + 
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size=16),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size=16),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16),
    strip.text.x = element_text(size=16),
    strip.text.y = element_text(size=16))

norm.viz
ggsave("./figures/norm_data_lineplot.svg", norm.viz)
