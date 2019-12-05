# Figure 3: Predicted biomass allocation

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

### Figure 3----
### 2.1.2 OLS (Unbiased)----
subti <- "OLS - unbiased"
dc <- melt(dc.ols$unbiased[, 1:5], id.vars = "xint")
dc$group <- rep(c("primary", "secondary"), each = nrow(dc)/2)

fig3_v1 <-
  ggplot(dc, aes(color = variable, y = value, x = xint, lty = group)) +
  geom_line(size = 3, alpha = al) +
  scale_y_continuous(breaks = breaks, limits = c(.10,.55),
                     labels = labels) +
  scale_color_manual(values = col, name = "", 
                     labels = c("Male", "Female", "Petals", "Sepals")) +
  scale_linetype_manual(values = c(1,5)) +
  scale_x_log10(breaks =  xla, labels = xla) +
  theme_classic(base_size = 20) +
  labs(x = "Flower biomass (g)", y = "Predicted allocation (%)") +
  guides(lty = FALSE) +
  theme(legend.position = c(0.15,.9))

fig3_v1 

## Add pie charts--------
# Create a basic bar
dc$xint <- as.factor(dc$xint)
dc3 <- dc %>% filter(xint %in% c(0.00001, 0.001, 1)) %>%  droplevels()

levels(dc3$xint) <- c("0.00001g", "0.001g", "1g")

g2 <-
  ggplot(dc3, aes(x="", y=value, fill = factor(variable, 
                                               levels = c("male", "petals","sepals", "female")[4:1]))) + 
  geom_bar(stat="identity", size = .25, color = "white", show.legend = F, alpha = .8) +
  coord_polar("y", start=0) + 
  #geom_text(aes(label = paste0(round(value*100, digits = 2), "%")),
  #          position = position_stack(vjust = 0.5)) + 
  scale_fill_manual(values = col[c(1,3,4,2)[4:1]]) + 
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_size = 20) + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        title =  element_text(), strip.background = element_blank()) +
  facet_grid(~xint);g2 

fig_3AB <-
  plot_grid(fig3_v1, g2, ncol = 1, 
            rel_heights = c(2,1), 
            rel_widths = c(2,2), 
            labels = c("A", "B"),
            label_size = 18, 
            label_fontface = "plain",
            align = "vh", axis = "l"); fig_3AB

### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/Fig_3_predicted_biomass.pdf",
       device = "pdf",
       plot = fig_3AB, width = 6, height = 8)
### END----