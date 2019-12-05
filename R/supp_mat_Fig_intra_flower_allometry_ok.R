# Supplementary material
# Figure: Intra-flower allometries

### Packages:-------------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(smatr)

options(scipen = 10000)
rm(list = ls())

source("R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
# dataset
d   <- read.csv("data/processed/data_flower_biomass_partition.csv")

### Intra-flower allometries----
### 1. Pair-wise between flower organs----------
### 1.1 Female ~ Male----
gynand <- sma(log10(gyn) ~ log10(and), data = d, slope.test = 1)
summary(gynand)

### 1.2 Female ~ Petals----
gynpet <- sma(log10(gyn) ~ log10(pet), data = d, slope.test = 1)
summary(gynpet)

### 1.3 Female ~ Sepals----
gynsep <- sma(log10(gyn) ~ log10(sep), data = d, slope.test = 1)
summary(gynsep)

### 1.4 Sepals ~ Petals----
seppet <- sma(log10(sep) ~ log10(pet), data = d, slope.test = 1)
summary(seppet)

### 1.5 Male ~ Petals----
andpet <- sma(log10(and) ~ log10(pet), data = d, slope.test = 1)
summary(andpet)

### 1.6 Male ~ Sepals----
andsep <- sma(log10(and) ~ log10(sep), data = d, slope.test = 1)
summary(andsep)

### Summarise tables------
stat_table1 <- 
  tab.rma(mod.list = list(gynand, 
                              gynpet,
                              gynsep, 
                              seppet,
                              andpet,
                              andsep),
              title = "phySMA regression", 
              model.names = c("Female ~ Male", 
                              "Female ~ Petals",
                              "Female ~ Sepals",
                              "Sepals ~ Petals",
                              "Male ~ Petals",
                              "Male ~ Sepals"))

### 2. Sexual ~ Flower biomass----------
### 2.1 Primary ~ Flower----
prim <- sma(log10(and + gyn) ~ log10(tot), data = d, slope.test = 1)
summary(prim)

### 2.2 Secondary ~ Flower----
seco <- sma(log10(pet + sep) ~ log10(tot), data = d, slope.test = 1)
summary(seco)

### 3. Primary sexual ~ secondary sexual organs----------
### 3.1 Female ~ Male----
cont <- sma(log10(and + gyn) ~ log10(pet + sep), data = d, slope.test = 1)
summary(cont)

### Summarise tables------
stat_table2 <- 
  tab.rma(mod.list = list(prim, seco, cont),
          title = "phySMA regression", 
          model.names = c("Primary organs ~ Flower", 
                          "Secondary organs ~ Flower",
                          "Primary organs ~ Secondary organs"))


# Figures----------------
# Pair-wise flower organs----
limx <- c(0.000005, 1)
intx <- c(0.000001 ,0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
limy <- c(0.000005, 1)
inty <- c(0.000001 ,0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
xa   <- 0.005 
ya   <- 0.0005
bsize <- 14
cols <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")

### 1.1 female ~ male----
g1 <-  plot_phylolm(y = "gyn", x = "and", data = d, estimates = stat_table1, 
                      subtitle = "Female ~ Male", 
                      laby = "Female biomass (g)", 
                      labx = "Male biomass (g)",
                      xa = xa,
                      ya = ya,
                      intx = intx, inty = inty,
                      limx = limx, limy = limy,
                      c.ponto = "steelblue", 
                      c.line = gray(.2), bsize = bsize)
g1

### 1.2 female ~ petals----
g2 <-  plot_phylolm(y = "gyn", x = "pet", data = d, estimates = stat_table1, 
                    subtitle = "Female ~ Petals", 
                    laby = "Female biomass (g)", 
                    labx = "Petals biomass (g)",
                    xa = xa,
                    ya = ya,
                    intx = intx, inty = inty,
                    limx = limx, limy = limy,
                    c.ponto = "steelblue", 
                    c.line = gray(.2), bsize = bsize)
g2

### 1.3 female ~ sepals----
g3 <-  plot_phylolm(y = "gyn", x = "sep", data = d, estimates = stat_table1, 
                    subtitle = "Female ~ Sepals", 
                    laby = "Female biomass (g)", 
                    labx = "Sepals biomass (g)",
                    xa = xa,
                    ya = ya,
                    intx = intx, inty = inty,
                    limx = limx, limy = limy,
                    c.ponto = "steelblue", 
                    c.line = gray(.2), bsize = bsize)
g3

### 1.4 sepals ~ petals----
g4 <-  plot_phylolm(y = "sep", x = "pet", data = d, estimates = stat_table1, 
                    subtitle = "Sepals ~ Petals", 
                    laby = "Sepals biomass (g)", 
                    labx = "Petals biomass (g)",
                    xa = xa,
                    ya = ya,
                    intx = intx, inty = inty,
                    limx = limx, limy = limy,
                    c.ponto = "steelblue", 
                    c.line = gray(.2), bsize = bsize)
g4

### 1.5 male ~ petals----
g5 <-  plot_phylolm(y = "and", x = "pet", data = d, estimates = stat_table1, 
                    subtitle = "Male ~ Petals", 
                    laby = "Male biomass (g)", 
                    labx = "Petals biomass (g)",
                    xa = xa,
                    ya = ya,
                    intx = intx, inty = inty,
                    limx = limx, limy = limy,
                    c.ponto = "steelblue", 
                    c.line = gray(.2), bsize = bsize)
g5

### 1.6 male ~ sepals----
g6 <-  plot_phylolm(y = "and", x = "sep", data = d, estimates = stat_table1, 
                    subtitle = "Male ~ Sepals", 
                    laby = "Male biomass (g)", 
                    labx = "Sepals biomass (g)",
                    xa = xa,
                    ya = ya,
                    intx = intx, inty = inty,
                    limx = limx, limy = limy,
                    c.ponto = "steelblue", 
                    c.line = gray(.2), bsize = bsize)
g6

gpair <- 
  plot_grid(g1, g2, g3, g4, g5, g6, labels = LETTERS[1:6],
          label_fontface = "plain", label_size = 18)
gpair

### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/supp/SFig_intra_flower_flower_organs.png",
       height = 7, width = 13, units = "in",
       plot = gpair)

# Primary versus secondary
limx <- c(0.00005, 1)
intx <- c(0.000001 ,0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
limy <- c(0.00001, 1)
inty <- c(0.000001 ,0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
xa   <- 0.05 
ya   <- 0.005
bsize <- 18
cols <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")

### 1.1 Primary ~ Flower----
d$pri <- d$and + d$gyn

g1 <-  plot_phylolm(y = "pri", x = "tot", data = d, estimates = stat_table2, 
                    subtitle = "Primary organs ~ Flower", 
                    laby = "Primary organs (g)", 
                    labx = "Flower biomass (g)",
                    xa = xa,
                    ya = ya,
                    intx = intx, inty = inty,
                    limx = limx, limy = limy,
                    c.ponto = "steelblue", 
                    c.line = gray(.2), bsize = bsize)
g1

### 1.2 Secondary ~ Flower----
d$sec <- d$pet + d$sep

g2 <-  plot_phylolm(y = "sec", x = "tot", data = d, estimates = stat_table2, 
                    subtitle = "Secondary organs ~ Flower", 
                    laby = "Secondary organs (g)", 
                    labx = "Flower biomass (g)",
                    xa = xa,
                    ya = ya,
                    intx = intx, inty = inty,
                    limx = limx, limy = limy,
                    c.ponto = "steelblue", 
                    c.line = gray(.2), bsize = bsize)
g2

### 1.3 primary ~ secondary----
g3 <-  plot_phylolm(y = "pri", x = "sec", data = d, estimates = stat_table2, 
                    subtitle = "Primary organs ~ Secondary organs", 
                    laby = "Primary organs (g)", 
                    labx = "Secondary organs (g)",
                    xa = xa,
                    ya = ya,
                    intx = intx, inty = inty,
                    limx = limx, limy = limy,
                    c.ponto = "steelblue", 
                    c.line = gray(.2), bsize = bsize)
g3

gcont <- 
  plot_grid(g1, g2, g3, labels = LETTERS[1:3], nrow = 3,
            label_fontface = "plain", label_size = 18)
gcont
ggsave("outputs/figures/supp/SFig_intra_flower_prim_vs_seco.png",
       height = 12, width = 6, units = "in",
       plot = gcont)
