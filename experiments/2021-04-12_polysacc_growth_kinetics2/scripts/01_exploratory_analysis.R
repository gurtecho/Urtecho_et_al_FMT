library(tidyverse)

rm(list=ls())

# Working directory for the experiment
setwd("~/Dropbox/projects/fmt/experiments/2021-04-12_polysacc_growth_kinetics2/")

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
  filter(polysacc == "Water") %>%
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

data

# convert sample type to factor which levels fixed in order I want to plot
data$sample_type = factor(data$sample_type, levels = c("Jax", "Envigo"))


# Order I want to plot polysacc in, to be used as levels
polysacc.order <-c("Water", "Arabinan", "Arabinogalactan", "A_Cellulose", "Cellulose",
                   "Dextran", "Mucin", "Pectic Galactan", "Pectin", "MCF")

# convert polysacc to factor, levels in order I want to plot
data$polysacc <- factor(data$polysacc, levels = polysacc.order)

data$replicate <- factor(data$replicate, levels=as.character(unique(data$replicate)))

# Plot raw data
raw.viz <- data  %>%
  ggplot() +
  geom_line(aes(x=time_mins, y=OD600, color=replicate)) +
  facet_grid(polysacc~sample_type) +
  xlab("Time (mins)") +
  labs(color="Replicate") +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size=16),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size=16),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16),
    strip.text.x = element_text(size=12),
    strip.text.y = element_text(size=12),
    strip.text.y.right = element_text(angle = 0))

raw.viz
ggsave("./figures/raw_data_lineplot.svg", raw.viz)
ggsave("./figures/raw_data_lineplot.pdf", raw.viz, width=8, height=10)
  

# Plot normalized data
norm.viz <- data  %>%
  filter(polysacc != "Water") %>%
  ggplot() +
  geom_line(aes(x=time_mins, y=norm_OD600, color=replicate)) +
  facet_grid(polysacc~sample_type) + 
  xlab("Time (mins)") +
  ylab("Normalized OD600") +
  labs(color="Replicate") +
  ylim(0, 4) + 
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size=12),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    strip.text.x = element_text(size=12),
    strip.text.y = element_text(size=12),
    strip.text.y.right = element_text(angle = 0))

norm.viz
ggsave("./figures/norm_data_lineplot.svg", norm.viz)
ggsave("./figures/norm_data_lineplot.pdf", norm.viz, width=8, height=10)


# Create averaged out data
avg.data <- data %>%
  group_by(polysacc, sample_type, time_mins) %>%
  summarise(OD600 = mean(OD600), bg= mean(bg))

# Compute OD600 normalized against water controls
avg.data$norm_OD600 <- avg.data$OD600 / avg.data$bg

# Compute the OD600 corrected against water controls
avg.data$corr_OD600 <- avg.data$OD600 - avg.data$bg
  

# Plot averaged out data
avg.viz <- avg.data  %>%
  ggplot() +
  geom_line(aes(x=time_mins, y=OD600, color=polysacc)) +
  facet_grid(polysacc~sample_type) +
  xlab("Time (mins)") +
  labs(color="Replicate") +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size=16),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size=16),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16),
    strip.text.x = element_text(size=12),
    strip.text.y = element_text(size=12),
    strip.text.y.right = element_text(angle = 0))

avg.viz
ggsave("./figures/avg_data_lineplot.svg", avg.viz)
ggsave("./figures/avg_data_lineplot.pdf", avg.viz, width=8, height=10)


# plot out averaged normalized data
# Plot normalized data
avg.norm.viz <- avg.data  %>%
  filter(polysacc != "Water") %>%
  ggplot() +
  geom_line(aes(x=time_mins, y=norm_OD600, color=polysacc)) +
  facet_grid(polysacc~sample_type) + 
  xlab("Time (mins)") +
  ylab("Normalized OD600") +
  labs(color="Replicate") +
  ylim(0, 4) + 
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size=12),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    strip.text.x = element_text(size=12),
    strip.text.y = element_text(size=12),
    strip.text.y.right = element_text(angle = 0))

avg.norm.viz
ggsave("./figures/norm_data_lineplot.svg", avg.norm.viz)
ggsave("./figures/norm_data_lineplot.pdf", avg.norm.viz, width=8, height=10)

avg.data

# Plot out avg data w/ Jax vs. Envigo together
grouped.avg.viz <- avg.data  %>%
  ggplot() +
  geom_line(aes(x=time_mins, y=OD600, color=sample_type)) +
  facet_grid(polysacc~.) +
  xlab("Time (mins)") +
  labs(color="Mouse Type") +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size=16),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size=16),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16),
    strip.text.x = element_text(size=12),
    strip.text.y = element_text(size=12),
    strip.text.y.right = element_text(angle = 0))

grouped.avg.viz
ggsave("./figures/grouped_avg_data_lineplot.svg", grouped.avg.viz)
ggsave("./figures/grouped_avg_data_lineplot.pdf", grouped.avg.viz, width=8, height=10)


