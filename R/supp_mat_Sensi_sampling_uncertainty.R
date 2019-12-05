# Supp. Mat - Sensitivity analysis (sampling uncertainty)
# Warning: Long time to run (~ 2 hours)

### Packages--------------------------------------------------------------------
library(sensiPhy)
library(phylolm)
library(tidyverse)
library(smatr)
library(phytools)
library(cowplot)

### Start-----------------------------------------------------------------------
rm(list = ls())
source("R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
p  <- read.tree("data/processed/working_study_phylo_tree.tre")
d  <- read.csv("data/processed/data_flower_biomass_partition.csv")
rownames(d) <- d$sp
cd <- match_dataphy(tot ~ 1, data = d, phy = p)
estimates <- read.csv("outputs/tables/supp/STable_model_stats_allometric_scalling.csv") 
  
### Sensitivity to Sampling Uncertaing----
### Parameters
times = 1000

### 1. Male ~ flower-------------------------------------------------------
### 1.1 SMA----
samp.andsma <- samp_sma(log10(and) ~ log10(tot), data = d, 
                    times = times, method = "SMA")
### 1.2 phySMA----
samp.andphysma <- samp_physma(cdata = cd, y = "and", x = "tot",  
                              method = "BM", times = times)
### 1.3 OLS----
samp.andols <- samp_sma(log10(and) ~ log10(tot), data = d, 
                     times = times, method = "OLS")
### 1.4 PGLS----
samp.andpgls <- samp_pgls(log10(and) ~ log10(tot), cdata = cd, 
                           times = times, method = "BM")

### 2. Female ~ flower------------------
### 2.1 SMA----
samp.gynsma <- samp_sma(log10(gyn) ~ log10(tot), data = d, 
                        times = times, method = "SMA")
### 2.2 phySMA----
samp.gynphysma <- samp_physma(cdata = cd, y = "gyn", x = "tot",  
                              method = "BM", times = times)
### 2.3 OLS----
samp.gynols <- samp_sma(log10(gyn) ~ log10(tot), data = d, 
                        times = times, method = "OLS")
### 2.4 PGLS----
samp.gynpgls <- samp_pgls(log10(gyn) ~ log10(tot), cdata = cd, 
                          times = times, method = "BM")

### 3. Petals ~ flower-------------------------------------------------------
### 3.1 SMA----
samp.petsma <- samp_sma(log10(pet) ~ log10(tot), data = d, 
                        times = times, method = "SMA")
### 3.2 phySMA----
samp.petphysma <- samp_physma(cdata = cd, y = "pet", x = "tot",  
                              method = "BM", times = times)
### 3.3 OLS----
samp.petols <- samp_sma(log10(pet) ~ log10(tot), data = d, 
                        times = times, method = "OLS")
### 3.4 PGLS----
samp.petpgls <- samp_pgls(log10(pet) ~ log10(tot), cdata = cd, 
                          times = times, method = "BM")

### 4. Sepals ~ flower-------------------------------------------------------
### 4.1 SMA----
samp.sepsma <- samp_sma(log10(sep) ~ log10(tot), data = d, 
                        times = times, method = "SMA")
### 4.2 phySMA----
samp.sepphysma <- samp_physma(cdata = cd, y = "sep", x = "tot",  
                              method = "BM", times = times)
### 4.3 OLS----
samp.sepols <- samp_sma(log10(sep) ~ log10(tot), data = d, 
                        times = times, method = "OLS")
### 4.4 PGLS----
samp.seppgls <- samp_pgls(log10(sep) ~ log10(tot), cdata = cd, 
                          times = times, method = "BM")

### Export results--------------------------------------------------------------
sensi_samp <- list(
  samp.andsma, samp.andphysma, samp.andols, samp.andpgls,
  samp.gynsma, samp.gynphysma, samp.gynols, samp.gynpgls,
  samp.petsma, samp.petphysma, samp.petols, samp.petpgls,
  samp.sepsma, samp.sepphysma, samp.sepols, samp.seppgls
  )
names(sensi_samp) <- c(
  "samp.andsma", "samp.andphysma", "samp.andols", "samp.andpgls",
  "samp.gynsma", "samp.gynphysma", "samp.gynols", "samp.gynpgls",
  "samp.petsma", "samp.petphysma", "samp.petols", "samp.petpgls",
  "samp.sepsma", "samp.sepphysma", "samp.sepols", "samp.seppgls"
)

saveRDS(object = sensi_samp, file = "outputs/temp/sensi_sampling_models.Rds")


# Plot Figures------------------------------------------------------------------
# Load cach
sensi_samp <- readRDS("outputs/temp/sensi_sampling_models.Rds")
estimates  <- read.csv("outputs/tables/supp/STable_model_stats_allometric_scalling.csv")

cols <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")
lims <- c(0.8, 1.3)
bre  <- seq(0.8, 1.3, .05)

# Plots - SMA---------------------------------------------------------------------------
# Male
ti <- "RMA regression"
st <- "Male ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[1]

g1 <- plot_sensi_samp(x = sensi_samp$samp.andsma, 
                      color = color, fill = fill, 
                      st = st) 
g1
# Female
st <- "Female ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[2]

g2 <- plot_sensi_samp(x = sensi_samp$samp.gynsma, 
                      color = color, fill = fill, 
                      st = st);g2

# Petals
st <- "Petals ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[3]

g3 <- plot_sensi_samp(x = sensi_samp$samp.petsma, 
                      color = color, fill = fill, 
                      st = st);g3
# Sepals
st <- "Sepals ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[4]

g4 <- plot_sensi_samp(x = sensi_samp$samp.sepsma, 
                      color = color, fill = fill, 
                      st = st);g4

# Plot grid
gglobal <- 
  plot_grid(plotlist = list(g1, g2, g3, g4), 
            labels = c("A", "B", "C", "D"),label_size = 18,label_fontface = "plain",
            ncol = 2, nrow = 2, align = "hv")
title <- ggdraw() + 
  draw_label("Sensitivity analysis - Sampling uncertainty (SMA)", size = 18)

gglobal <- 
  plot_grid(title, gglobal, ncol = 1, rel_heights = c(0.1, 1))

### Save plot
ggsave("outputs/figures/supp/SFig_sensi_sampling_sma_regressions.png",
       height = 9, width = 10, units = "in",
       plot = gglobal)

# Plots - phySMA---------------------------------------------------------------------------
lims <- c(0.8, 1.4)
bre  <- seq(0.8, 1.4, .05)

# Male
ti <- "RMA.phy regression"
st <- "Male ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[1]

g1 <- plot_sensi_samp(x = sensi_samp$samp.andphysma, 
                      color = color, fill = fill, 
                      st = st);g1

# Female
st <- "Female ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[2]

g2 <- plot_sensi_samp(x = sensi_samp$samp.gynphysma, 
                      color = color, fill = fill, 
                      st = st);g2

# Petals
st <- "Petals ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[3]

g3 <- plot_sensi_samp(x = sensi_samp$samp.petphysma, 
                      color = color, fill = fill, 
                      st = st);g3


# Sepals
st <- "Sepals ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[4]

g4 <- plot_sensi_samp(x = sensi_samp$samp.sepphysma, 
                      color = color, fill = fill, 
                      st = st);g4

# Plot grid
gglobal <- 
  plot_grid(plotlist = list(g1, g2, g3, g4), 
            labels = c("A", "B", "C", "D"),label_size = 18,label_fontface = "plain",
            ncol = 2, nrow = 2, align = "hv")
title <- ggdraw() + 
  draw_label("Sensitivity analysis - Sampling uncertainty (phySMA)", size = 18)

gglobal <- 
  plot_grid(title, gglobal, ncol = 1, rel_heights = c(0.1, 1))

### Save plot
ggsave("outputs/figures/supp/SFig_sensi_sampling_physma_regressions.png",
       height = 9, width = 10, units = "in",
       plot = gglobal)

# Plots - OLS---------------------------------------------------------------------------
lims <- c(0.7, 1.25)
bre  <- seq(0.7, 1.25, .05)

# Male
ti <- "OLS regression"
st <- "Male ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[1]

g1 <- plot_sensi_samp(x = sensi_samp$samp.andols, 
                      color = color, fill = fill, 
                      st = st)

# Female
st <- "Female ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[2]

g2 <- plot_sensi_samp(x = sensi_samp$samp.gynols, 
                      color = color, fill = fill, 
                      st = st);g2

# Petals
st <- "Petals ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[3]

g3 <- plot_sensi_samp(x = sensi_samp$samp.petols, 
                      color = color, fill = fill, 
                      st = st);g3
# Sepals
st <- "Sepals ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[4]

g4 <- plot_sensi_samp(x = sensi_samp$samp.sepols, 
                      color = color, fill = fill, 
                      st = st);g4

# Plot grid
gglobal <- 
  plot_grid(plotlist = list(g1, g2, g3, g4), 
            labels = c("A", "B", "C", "D"),label_size = 18,label_fontface = "plain",
            ncol = 2, nrow = 2, align = "hv")
title <- ggdraw() + 
  draw_label("Sensitivity analysis - Sampling uncertainty (OLS)", size = 18)

gglobal <- 
  plot_grid(title, gglobal, ncol = 1, rel_heights = c(0.1, 1))

### Save plot
ggsave("outputs/figures/supp/SFig_sensi_sampling_ols_regressions.png",
       height = 9, width = 10, units = "in",
       plot = gglobal)

# Plots - PGLS---------------------------------------------------------------------------
# Male
ti <- "PGLS regression"
st <- "Male ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[1]

g1 <- plot_sensi_samp(x = sensi_samp$samp.andpgls, 
                      color = color, fill = fill, 
                      st = st)

# Female
st <- "Female ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[2]

g2 <- plot_sensi_samp(x = sensi_samp$samp.gynpgls, 
                      color = color, fill = fill, 
                      st = st);g2

# Petals
st <- "Petals ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[3]

g3 <- plot_sensi_samp(x = sensi_samp$samp.petpgls, 
                      color = color, fill = fill, 
                      st = st);g3
# Sepals
st <- "Sepals ~ Flower"
es <- estimates %>% filter(Title == ti & Subtitle == st)
color <- fill <- cols[4]

g4 <- plot_sensi_samp(x = sensi_samp$samp.seppgls, 
                      color = color, fill = fill, 
                      st = st);g4

# Plot grid
gglobal <- 
  plot_grid(plotlist = list(g1, g2, g3, g4), 
            labels = c("A", "B", "C", "D"),label_size = 18,label_fontface = "plain",
            ncol = 2, nrow = 2, align = "hv")
title <- ggdraw() + 
  draw_label("Sensitivity analysis - Sampling uncertainty (PGLS)", size = 18)

gglobal <- 
  plot_grid(title, gglobal, ncol = 1, rel_heights = c(0.1, 1))

### Save plot
ggsave("outputs/figures/supp/SFig_sensi_sampling_pgls_regressions.png",
       height = 9, width = 10, units = "in",
       plot = gglobal)
### END-----