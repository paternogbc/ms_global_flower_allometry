# Supplementary material
# Figure: Allometric scalling with axis independence (part ~ (total - part))

### Packages:-------------------------------------------------------------------
library(tidyverse)
library(cowplot)

options(scipen = 10000)
rm(list = ls())

source("R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
# dataset
d   <- read.csv("data/processed/data_flower_biomass_partition.csv")

# SMA------
est  <- read_csv("outputs/tables/supp/STable_sensi_axis_model_stats_allometric_scalling.csv")
est  <- filter(est, Title == "RMA regression (y ~ [x - Y])")

### 1. Scaling allometries (pgls)--------------------------------------------
### Define plot axis limits and breaks
limx <- c(0.00005, 1)
intx <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
limy <- c(0.000005, 1)
inty <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
xa   <- 0.05 
ya   <- 0.0005

cols <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")

### 1.1 male ~ flower----
dc <- d 
dc$tot <- dc$tot - dc$and

gand <-  plot_phylolm(y = "and", x = "tot", data = dc, estimates = est, 
                      subtitle = "Male ~ Flower", 
                      labx = "(Flower - Male) biomass (g)", 
                      laby = "Male biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[1], 
                      c.line = gray(.2))
gand

### 1.2 female ~ flower----
dc <- d 
dc$tot <- dc$tot - dc$gyn

ggyn <-  plot_phylolm(y = "gyn", x = "tot", data = dc, estimates = est, 
                      subtitle = "Female ~ Flower", 
                      labx = "(Flower - Female) biomass (g)", 
                      laby = "Female biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[2], 
                      c.line = gray(.2))
ggyn

### 1.3 petals ~ flower----
dc <- d 
dc$tot <- dc$tot - dc$pet

gpet <-  plot_phylolm(y = "pet", x = "tot", data = dc, estimates = est, 
                      subtitle = "Petals ~ Flower", 
                      labx = "(Flower - Petals) biomass (g)", 
                      laby = "Petals biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[3], 
                      c.line = gray(.2))
gpet

### 1.4 sepals ~ flower----
dc <- d 
dc$tot <- dc$tot - dc$sep

gsep <-  plot_phylolm(y = "sep", x = "tot", data = dc, estimates = est, 
                      subtitle = "Sepals ~ Flower", 
                      labx = "(Flower - Sepals) biomass (g)", 
                      laby = "Sepals biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[4], 
                      c.line = gray(.2))
gsep

## Add flower iiluestrations----------------------------------------------------
gand2 <-
  ggdraw() +
  draw_plot(gand) +
  draw_image(image = "images/flower_ilustrations/male.png", 
             x = -.15, 
             y = .32, 
             scale = .2)
ggyn2 <-
  ggdraw() +
  draw_plot(ggyn) +
  draw_image(image = "images/flower_ilustrations/female.png", 
             x = -.15, 
             y = .32, 
             scale = .2)
gpet2 <-
  ggdraw() +
  draw_plot(gpet) +
  draw_image(image = "images/flower_ilustrations/petals.png", 
             x = -.15, 
             y = .32, 
             scale = .2)
gsep2 <-
  ggdraw() +
  draw_plot(gsep) +
  draw_image(image = "images/flower_ilustrations/sepals.png", 
             x = -.15, 
             y = .32, 
             scale = .2)

### Figure ----
gglobal <- 
  plot_grid(plotlist = list(gand2, ggyn2, gpet2, gsep2), 
            labels = c("A", "B", "C", "D"),label_size = 18,label_fontface = "plain",
            ncol = 2, nrow = 2);gglobal

title <- ggdraw() + 
  draw_label("SMA allometries controlling for axis non-independence (Y ~ (X - Y))", size = 18)

gglobal <- 
  plot_grid(title, gglobal, ncol = 1, rel_heights = c(0.1, 1))
gglobal

### Save plot
ggsave("outputs/figures/supp/SFig_sensi_axis_global_allometry_sma.png",
       height = 8.27, width = 11.69, units = "in",
       plot = gglobal)

# phySMA------
est  <- read_csv("outputs/tables/supp/STable_sensi_axis_model_stats_allometric_scalling.csv")
est  <- filter(est, Title == "RMA.phy regression (y ~ [x - Y])")

### 1.1 male ~ flower----
dc <- d 
dc$tot <- dc$tot - dc$and

gand <-  plot_phylolm(y = "and", x = "tot", data = dc, estimates = est, 
                      subtitle = "Male ~ Flower", 
                      labx = "(Flower - Male) biomass (g)", 
                      laby = "Male biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[1], 
                      c.line = gray(.2))
gand

### 1.2 female ~ flower----
dc <- d 
dc$tot <- dc$tot - dc$gyn

ggyn <-  plot_phylolm(y = "gyn", x = "tot", data = dc, estimates = est, 
                      subtitle = "Female ~ Flower", 
                      labx = "(Flower - Female) biomass (g)", 
                      laby = "Female biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[2], 
                      c.line = gray(.2))
ggyn

### 1.3 petals ~ flower----
dc <- d 
dc$tot <- dc$tot - dc$pet

gpet <-  plot_phylolm(y = "pet", x = "tot", data = dc, estimates = est, 
                      subtitle = "Petals ~ Flower", 
                      labx = "(Flower - Petals) biomass (g)", 
                      laby = "Petals biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[3], 
                      c.line = gray(.2))
gpet

### 1.4 sepals ~ flower----
dc <- d 
dc$tot <- dc$tot - dc$sep

gsep <-  plot_phylolm(y = "sep", x = "tot", data = dc, estimates = est, 
                      subtitle = "Sepals ~ Flower", 
                      labx = "(Flower - Sepals) biomass (g)", 
                      laby = "Sepals biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = cols[4], 
                      c.line = gray(.2))
gsep

## Add flower iiluestrations----------------------------------------------------
gand2 <-
  ggdraw() +
  draw_plot(gand) +
  draw_image(image = "images/flower_ilustrations/male.png", 
             x = -.15, 
             y = .32, 
             scale = .2)
ggyn2 <-
  ggdraw() +
  draw_plot(ggyn) +
  draw_image(image = "images/flower_ilustrations/female.png", 
             x = -.15, 
             y = .32, 
             scale = .2)
gpet2 <-
  ggdraw() +
  draw_plot(gpet) +
  draw_image(image = "images/flower_ilustrations/petals.png", 
             x = -.15, 
             y = .32, 
             scale = .2)
gsep2 <-
  ggdraw() +
  draw_plot(gsep) +
  draw_image(image = "images/flower_ilustrations/sepals.png", 
             x = -.15, 
             y = .32, 
             scale = .2)

### Figure ----
gglobal <- 
  plot_grid(plotlist = list(gand2, ggyn2, gpet2, gsep2), 
            labels = c("A", "B", "C", "D"),label_size = 18,
            label_fontface = "plain",
            ncol = 2, nrow = 2);gglobal

title <- ggdraw() + 
  draw_label("phySMA allometries controlling for axis non-independence (Y ~ (X - Y))", size = 18)

gglobal <- 
  plot_grid(title, gglobal, ncol = 1, rel_heights = c(0.1, 1))

### Save plot
ggsave("outputs/figures/supp/SFig_sensi_axis_global_allometry_physma.png",
       height = 8.27, width = 11.69, units = "in",
       plot = gglobal)
### END-----