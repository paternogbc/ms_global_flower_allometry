# Supplementary material
# Figure: Predicted biomass for all models

### Start-----------------------------------------------------------------------
rm(list = ls())
source("R/ZZZ_functions.R")

### Packages:-------------------------------------------------------------------
library(sensiPhy)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(cowplot)
options(scipen = 100)

### Load Data:
est  <- read_csv("outputs/tables/supp/STable_model_stats_allometric_scalling.csv")
d    <- read_csv("data/processed/data_flower_biomass_partition.csv")

### Backtransform intercepts-----
est <- 
  est %>%
  mutate(int = 10^int,
         int.lo = 10^int.lo,
         int.hi = 10^int.hi)

### 1. Predicted biomass allocation----
e.ols <- filter(est, Title == "OLS regression")
e.pgls <- filter(est, Title == "PGLS regression")
e.rma <- filter(est, Title == "RMA regression")

### Predic bias and unbiased flower allocation----
dc.ols  <- pred.flower(estimates = e.ols)
dc.pgls <- pred.flower(estimates = e.pgls)
dc.rma  <- pred.flower(estimates = e.rma)

### 2. Plot predicted biomass allocation----
col    <- c("#D55E00","#56B4E9", "#CC79A7", "#009E73")
xla    <- c(0.00001, 0.0001, 0.001, 0.01, 0.1,1)
al     <- .7
lsize  <- 2
bs     <- 14
breaks <- seq(0.05,.7,.05)
limits <- c(0.05, .7)
labels <- seq(0.05,.7,.05)*100

### 2.1.1 OLS (Biased)----
subti <- "OLS - biased"
dc <- melt(dc.ols$biased, id.vars = "xint")

gols1 <-
  ggplot(dc, aes(color = variable, y = value, x = xint)) +
  geom_line(size = lsize, alpha = al, show.legend = T) +
  scale_y_continuous(breaks = breaks, limits = limits,
                     labels = labels) +
  theme_classic(base_size = bs) +
  scale_color_manual(values = col, name = "") +
  scale_x_log10(breaks =  xla, labels = xla) +
  theme(legend.position = c(.2,.9)) +
  labs(x = "Flower biomass (g)", y = "Predicted allocation (%)",
       subtitle = subti)
gols1

### 2.1.2 OLS (Unbiased)----
subti <- "OLS - unbiased"
dc <- melt(dc.ols$unbiased[, 1:5], id.vars = "xint")

gols2 <-
  ggplot(dc, aes(color = variable, y = value, x = xint)) +
  geom_line(size = lsize, alpha = al, show.legend = F) +
  scale_y_continuous(breaks = breaks, limits = limits,
                     labels = labels) +
  theme_classic(base_size = bs) +
  scale_color_manual(values = col, name = "") +
  scale_x_log10(breaks =  xla, labels = xla) +
  labs(x = "Flower biomass (g)", y = "Predicted allocation (%)",
       subtitle = subti)
gols2


### 2.2.1 PGLS (Biased)----
subti <- "PGLS - biased"
dc <- melt(dc.pgls$biased, id.vars = "xint")

gpgls1 <-
  ggplot(dc, aes(color = variable, y = value, x = xint)) +
  geom_line(size = lsize, alpha = al, show.legend = F) +
  scale_y_continuous(breaks = breaks, limits = limits,
                     labels = labels) +
  theme_classic(base_size = bs) +
  scale_color_manual(values = col, name = "") +
  scale_x_log10(breaks =  xla, labels = xla) +
  theme(legend.position = c(.2,.9)) +
  labs(x = "Flower biomass (g)", y = "Predicted allocation (%)",
       subtitle = subti)
gpgls1

### 2.2.2 PGLS (Unbiased)----
subti <- "PGLS - unbiased"
dc <- melt(dc.pgls$unbiased[, 1:5], id.vars = "xint")

gpgls2 <-
  ggplot(dc, aes(color = variable, y = value, x = xint)) +
  geom_line(size = lsize, alpha = al, show.legend = F) +
  scale_y_continuous(breaks = breaks, limits = limits,
                     labels = labels) +
  theme_classic(base_size = bs) +
  scale_color_manual(values = col, name = "") +
  scale_x_log10(breaks =  xla, labels = xla) +
  labs(x = "Flower biomass (g)", y = "Predicted allocation (%)",
       subtitle = subti)
gpgls2


### All plots----
gall <- 
  plot_grid(gols1, gpgls1, gols2, gpgls2,
            labels = LETTERS[1:4], 
            label_fontface = "plain", label_size = 18)
gall

### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/supp/SFig_predicted_biomass.png",
       height = 8.27, width = 9, units = "in",
       plot = gall)

### 3. Bias Plot------ 
# 3.1.1 OLS bias----
dc <- melt(dc.ols$biased[, 1:5], id.vars = "xint")
g1 <- 
  ggplot(dc, aes(y = value, x = xint, fill = variable)) +
  geom_bar(stat = "identity", alpha = al, show.legend = T) +
  scale_y_continuous(breaks = seq(0,1.2,.2), limits = c(0, 1.2), 
                     labels = paste(seq(0,1.2,.2)*100, "%")) +
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_fill_manual(values = col, name = "") +
  scale_x_log10(breaks =  c(10^seq(-5,0)), labels = xla) +
  theme_classic(base_size = bs) +
  theme(legend.position = c(.5,.95),
        legend.text = element_text(size = 10),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(x = .5, units = "cm"),
        legend.key.width  = unit(x = .5, units = "cm")) +
  labs(x = "Flower biomass", y = "Predicted allocation (%)", 
       subtitle = "OLS - biased")
g1

# 3.1.2 OLS unbiased----
dc <- melt(dc.ols$unbiased[, 1:5], id.vars = "xint")
g2 <- 
  ggplot(dc, aes(y = value, x = xint, fill = variable)) +
  geom_bar(stat = "identity", alpha = al, show.legend = F) +
  scale_y_continuous(breaks = seq(0,1.25,.25),
                     labels = paste(seq(0,1.25,.25)*100, "%")) +
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_fill_manual(values = col) +
  scale_x_log10(breaks =  c(10^seq(-5,0)), labels = xla) +
  theme_classic(base_size = bs) +
  labs(x = "Flower biomass", y = "Predicted allocation (%)",
       subtitle = "OLS - unbiased")
g2

# 3.2.1 PGLS bias----
dc <- melt(dc.pgls$biased[, 1:5], id.vars = "xint")
g3 <- 
  ggplot(dc, aes(y = value, x = xint, fill = variable)) +
  geom_bar(stat = "identity", alpha = al, show.legend = F) +
  scale_y_continuous(breaks = seq(0,1.2,.2), limits = c(0, 1.2), 
                     labels = paste(seq(0,1.2,.2)*100, "%")) +
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_fill_manual(values = col) +
  scale_x_log10(breaks =  c(10^seq(-5,0)), labels = xla) +
  theme_classic(base_size = bs) +
  labs(x = "Flower biomass", y = "Predicted allocation (%)",
       subtitle = "PGLS - biased")
g3

# 3.2.2 PGLS unbiased----
dc <- melt(dc.pgls$unbiased[, 1:5], id.vars = "xint")
g4 <- 
  ggplot(dc, aes(y = value, x = xint, fill = variable)) +
  geom_bar(stat = "identity", alpha = al, show.legend = F) +
  scale_y_continuous(breaks = seq(0,1.25,.25),
                     labels = paste(seq(0,1.25,.25)*100, "%")) +
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_fill_manual(values = col) +
  scale_x_log10(breaks =  c(10^seq(-5,0)), labels = xla) +
  theme_classic(base_size = bs) +
  labs(x = "Flower biomass", y = "Predicted allocation (%)",
       subtitle = "PGLS - unbiased")
g4

### All plots----
gall <- plot_grid(g1, g3, g2, g4, labels = LETTERS[1:4], 
                  label_fontface = "plain", label_size = 18)
gall

### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/supp/SFig_predicted_bias.png",
       height = 8.27, width = 9, units = "in",
       plot = gall)
### END----