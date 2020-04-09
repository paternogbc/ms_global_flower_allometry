# Figure 2: Global allometric regressions (flower components ~ flower biomass)

### Packages:-------------------------------------------------------------------
library(tidyverse)
library(cowplot)

options(scipen = 10000)
rm(list = ls())

source("R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
# dataset
d   <- read.csv("data/processed/data_flower_biomass_partition.csv")

# SMA model estimates---
est  <- read_csv("outputs/tables/supp/STable_model_stats_allometric_scalling.csv")
est  <- filter(est, Title == "RMA regression")

### 1. Scaling allometries (pgls)--------------------------------------------
### Define plot axis limits and breaks
limx  <- c(0.00005, 1)
intx  <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
limy  <- c(0.000005, 1)
inty  <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
xa    <- 0.028 
ya    <- 0.00028
bsize <- 12
anota.size <- 3.2
axis_size  <- 14

cols <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")

### 1.1 male ~ flower----
gand <-  plot_phylolm(y = "and", x = "tot", data = d, estimates = est, 
                      subtitle = "Male ~ Flower", 
                      labx = "Flower biomass (g)", 
                      laby = "Male biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[1], 
                      c.line = gray(.2),
                      bsize = bsize,
                      anota.size = anota.size,
                      axis_size = axis_size)
gand

### 1.2 female ~ flower----
ggyn <-  plot_phylolm(y = "gyn", x = "tot", data = d, estimates = est, 
                      subtitle = "Female ~ Flower", 
                      labx = "Flower biomass (g)", 
                      laby = "Female biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[2], 
                      c.line = gray(.2),
                      bsize = bsize,
                      anota.size = anota.size,
                      axis_size = axis_size)
ggyn

### 1.3 petals ~ flower----
gpet <-  plot_phylolm(y = "pet", x = "tot", data = d, estimates = est, 
                      subtitle = "Petals ~ Flower", 
                      labx = "Flower biomass (g)", 
                      laby = "Petals biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[3], 
                      c.line = gray(.2),
                      bsize = bsize,
                      anota.size = anota.size,
                      axis_size = axis_size)
gpet

### 1.4 sepals ~ flower----
gsep <-  plot_phylolm(y = "sep", x = "tot", data = d, estimates = est, 
                      subtitle = "Sepals ~ Flower", 
                      labx = "Flower biomass (g)", 
                      laby = "Sepals biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[4], 
                      c.line = gray(.2),
                      bsize = bsize,
                      anota.size = anota.size,
                      axis_size = axis_size)
gsep

## Add flower iiluestrations----------------------------------------------------
gand2 <-
  ggdraw() +
  draw_plot(gand) +
  draw_image(image = "images/flower_ilustrations/male.png", 
             x = -.1, 
             y = .34, 
             scale = .23)
ggyn2 <-
  ggdraw() +
  draw_plot(ggyn) +
  draw_image(image = "images/flower_ilustrations/female.png", 
             x = -.1, 
             y = .34, 
             scale = .23)
gpet2 <-
  ggdraw() +
  draw_plot(gpet) +
  draw_image(image = "images/flower_ilustrations/petals.png", 
             x = -.1, 
             y = .34, 
             scale = .23)
gsep2 <-
  ggdraw() +
  draw_plot(gsep) +
  draw_image(image = "images/flower_ilustrations/sepals.png", 
             x = -.1, 
             y = .34, 
             scale = .23)
### Figure ----
gglobal <- 
  plot_grid(plotlist = list(gand2, ggyn2, gpet2, gsep2), 
            labels = c("A", "B", "C", "D"),label_size = 18,label_fontface = "plain",
            ncol = 2, nrow = 2);gglobal


### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/Fig_2_global_allometry_sma.pdf",
       device = "pdf",
       height = 15, width = 17.8, units = "cm",
       plot = gglobal)
### END-----