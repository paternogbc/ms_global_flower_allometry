# Supplementary material
# Figure: Taxonomic Coverage (Top orders and families)

# Packages ----------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(reshape2)
library(cowplot)
source(file = "R/ZZZ_functions.R")

# Load data ---------------------------------------------------------------
tax <- read.csv("outputs/tables/supp/STable_plant_names_taxonomy.csv")
tax <- select(tax, order, family)

### Distribution of species across Orders
dc <- count(tax, order)

g1 <- ggplot(dc, aes(y = n, x = reorder(order, -n))) +
  geom_col(fill = "lightgreen") +
  tema(base_size = 20) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, size = 14, hjust = 1)) +
  xlab("Order") + ylab("Number of species"); g1

### Distribution of species across Families
dc <- count(tax, family)
  
g2 <- ggplot(dc, aes(y = n, x = reorder(family, -n))) +
  geom_col(fill = "lightgreen") +
  tema(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Family") + ylab("Number of species"); g2


### Export graphs to A4 landscape format (pdf):
ggsave(plot = g1, filename = "outputs/figures/supp/SFig_taxonomic_orders.png",
       height = 5, width = 11.69, units = "in")
ggsave(plot = g2, filename = "outputs/figures/supp/SFig_taxonomic_families.png",
       height = 4, width = 11.69, units = "in")
