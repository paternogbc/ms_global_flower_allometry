# Supp. Mat
# Table with outlying species in terms of divergence from the fitted lines (SMA)

# Packages----------------------------------------------------------------------
library(tidyverse)
library(ggrepel)

# Start-------------------------------------------------------------------------
rm(list = ls())
source("R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
d <- read.csv("data/processed/data_flower_biomass_partition.csv")

# SMA model fits
m <- readRDS(file = "outputs/temp/fitted_models_allometric_scaling.Rds")

# Outlying species from fitted lines--------------------------------------------

# 1. Male ~ FLower biomass---------------
cols <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")

res  <- resid(m$rmaand)
sres <- (mean(res) - res)/sd(res)
spand <- m$rmaand$data[which(abs(sres) > 3),] %>% rownames()

# Check if reming these species alterns the main SMA regression
summary(sma(log10(and) ~ log10(tot), data = filter(d, !sp %in% spand), slope.test = 1))

g1 <-
  ggplot(d, aes(y = log10(and), x = log10(tot))) +
  geom_point(color = "black", shape = 21, fill = cols[1], alpha = .5) +
  geom_point(data = filter(d, sp %in% spand), color = "red", size = 3) +
  geom_label_repel(data = filter(d, sp %in% spand), aes(label = sp), 
                  color = "red") +
  tema(base_size = 16) +
  labs(x = "Flower biomass (g)[log]",
       y = "Male biomass (g)[log]",
       subtitle = "")
g1

# 2. Female ~ FLower biomass---------------
res  <- resid(m$rmagyn)
sres <- (mean(res) - res)/sd(res)
spgyn <- m$rmagyn$data[which(abs(sres) > 3),] %>% rownames()

# Check if reming these species alterns the main SMA regression
summary(sma(log10(gyn) ~ log10(tot), data = filter(d, !sp %in% spgyn), slope.test = 1))

g2 <-
  ggplot(d, aes(y = log10(gyn), x = log10(tot))) +
  geom_point(color = "black", shape = 21, fill = cols[2], alpha = .5) +
  geom_point(data = filter(d, sp %in% spgyn), color = "red", size = 3) +
  geom_label_repel(data = filter(d, sp %in% spgyn), aes(label = sp), 
                   color = "red") +
  tema(base_size = 16) +
  labs(x = "Flower biomass (g)[log]",
       y = "Female biomass (g)[log]")
g2

# 3. Petals ~ FLower biomass---------------
res  <- resid(m$rmapet)
sres <- (mean(res) - res)/sd(res)
sppet <- m$rmapet$data[which(abs(sres) > 3),] %>% rownames()

# Check if reming these species alterns the main SMA regression
summary(sma(log10(pet) ~ log10(tot), data = filter(d, !sp %in% sppet), slope.test = 1))

g3 <-
  ggplot(d, aes(y = log10(pet), x = log10(tot))) +
  geom_point(color = "black", shape = 21, fill = cols[3], alpha = .5) +
  geom_point(data = filter(d, sp %in% sppet), color = "red", size = 3) +
  geom_label_repel(data = filter(d, sp %in% sppet), aes(label = sp), 
                   color = "red") +
  tema(base_size = 16) +
  labs(x = "Flower biomass (g)[log]",
       y = "Petals biomass (g)[log]")
g3

# 4. Sepals ~ FLower biomass---------------
res  <- resid(m$rmasep)
sres <- (mean(res) - res)/sd(res)
spsep <- m$rmasep$data[which(abs(sres) > 3),] %>% rownames()

# Check if reming these species alterns the main SMA regression
summary(sma(log10(sep) ~ log10(tot), data = filter(d, !sp %in% spsep), slope.test = 1))

g4 <-
  ggplot(d, aes(y = log10(sep), x = log10(tot))) +
  geom_point(color = "black", shape = 21, fill = cols[4], alpha = .5) +
  geom_point(data = filter(d, sp %in% spsep), color = "red", size = 3) +
  geom_label_repel(data = filter(d, sp %in% spsep), aes(label = sp), 
                   color = "red", force = 2) +
  tema(base_size = 16) +
  labs(x = "Flower biomass (g)[log]",
       y = "sepals biomass (g)[log]")

# Plot grid
gglobal <- 
  plot_grid(plotlist = list(g1, g2, g3, g4), 
            labels = c("A", "B", "C", "D"),label_size = 18,label_fontface = "plain",
            ncol = 2, nrow = 2, align = "hv")
title <- ggdraw() + 
  draw_label("Outlying species from SMA fitted lines (std residuals > 3)", size = 18)

gglobal <- 
  plot_grid(title, gglobal, ncol = 1, rel_heights = c(0.1, 1))

### Save plot----------------
ggsave("outputs/figures/supp/SFig_Sensi_outlying_species_sma_regressions.png",
       height = 9, width = 10, units = "in",
       plot = gglobal)
# END--------------