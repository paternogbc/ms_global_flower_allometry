# Supp. Mat - Sensitivity analysis (Phylogenetic uncertainty)

### Packages:-------------------------------------------------------------------
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
p   <- read.tree("data/processed/300_study_phylo_trees.tre")
d   <- read.csv("data/processed/data_flower_biomass_partition.csv")
estimates <- read.csv("outputs/tables/supp/STable_model_stats_allometric_scalling.csv") 

# Match with Phylogeny----------------------------------------------------------
rownames(d) <- d$sp
d <- match_dataphy(tot ~ 1, data = d, phy = p)

### 1. Male ~ Flower----
### phySMA----
and_sma <- tree_physma(y = "and", x = "tot", data = d$data, trees = d$phy, method = "BM")

### PGLS----
and_pgls <- tree_phylm(log10(and) ~ log10(tot), data = d$data, phy = d$phy,
                       n.tree = 300, model = "BM")
### 2. Female ~ Flower----
### phySMA----
gyn_sma <- tree_physma(y = "gyn", x = "tot", data = d$data, trees = d$phy, 
                       method = "BM")

### PGLS----
gyn_pgls <- tree_phylm(log10(gyn) ~ log10(tot), data = d$data, phy = d$phy,
                       n.tree = 300, model = "BM")

### 3. Petals ~ Flower----
### phySMA----
pet_sma <- tree_physma(y = "pet", x = "tot", data = d$data, trees = d$phy, method = "BM")

### PGLS----
pet_pgls <- tree_phylm(log10(pet) ~ log10(tot), data = d$data, phy = d$phy,
                       n.tree = 300, model = "BM")

### 4. Sepals ~ Flower----
### phySMA----
sep_sma <- tree_physma(y = "sep", x = "tot", data = d$data, trees = d$phy, method = "BM")

### PGLS----
sep_pgls <- tree_phylm(log10(sep) ~ log10(tot), data = d$data, phy = d$phy,
                       n.tree = 300, model = "BM")



# Save sensitivity cache--------------------------------------------------------
sensi <- list("and_physma" = and_sma,
              "and_pgls" = and_pgls$sensi.estimates,
              "gyn_physma" = gyn_sma,
              "gyn_pgls" = gyn_pgls$sensi.estimates,
              "pet_physma" = pet_sma,
              "pet_pgls" = pet_pgls$sensi.estimates,
              "sep_physma" = sep_sma,
              "sep_pgls" = sep_pgls$sensi.estimates)

saveRDS(object = sensi, file = "outputs/temp/sensi_phylogenetic_uncertainty.Rds")

# Plot Figures------------------------------------------------------------------
sensi <- readRDS("outputs/temp/sensi_phylogenetic_uncertainty.Rds")
cols <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")
lims <- c(0.8, 1.3)
bre1  <- seq(0.8, 1.3, .05)

# phySMA-----
est_physma <- data.frame(
  estimate = c(sensi$and_physma$estimate, 
               sensi$gyn_physma$estimate,
               sensi$pet_physma$estimate,
               sensi$sep_physma$estimate),
  organ = rep(c("and", "gyn", "pet", "sep"), each = nrow(sensi$and_physma)))

g1 <-
  ggplot(est_physma %>% filter(organ %in% c("and", "gyn")), aes(x = estimate, fill = organ)) +
  geom_histogram(position = position_dodge(), 
                 bins = 20, show.legend = T,
                 alpha = .8) +
  scale_fill_manual(values = cols[1:2], name = "", labels = c("Male", "Female")) +
  scale_x_continuous(breaks = seq(0.9, 1.3, 0.05)) +
  labs(y = "Frequency", x = "Estimated slopes",
       subtitle = "phySMA | (Male vs Female)") +
  tema(base_size = 20) +
  theme(legend.position = c(.88,.9),
        legend.background = element_blank());g1
  
g2 <-
  ggplot(est_physma %>% filter(organ %in% c("pet", "sep")),
       aes(x = estimate, fill = organ)) +
  geom_histogram(position = position_dodge(), 
                 bins = 20, show.legend = T,
                 alpha = .8) +
  scale_fill_manual(values = cols[3:4],
                    name = "", labels = c("Petals", "Sepals")) +
  scale_x_continuous(breaks = seq(0.8, 1.5, 0.1)) +
  labs(y = "Frequency", x = "Estimated slopes",
       subtitle = "phySMA | (Petals vs Sepals)") +
  tema(base_size = 20) +
  theme(legend.position = c(.88,.9),
        legend.background = element_blank());g2

gsma <- plot_grid(g1, g2,labels = LETTERS[1:2], 
                  label_size = 20, label_fontface = "plain"); gsma

title <- ggdraw() + 
  draw_label("Sensitivity analysis - Phylogenetic uncertainty (phySMA)", size = 18)

gglobal <- 
  plot_grid(title, gsma, ncol = 1, rel_heights = c(0.1, 1));gglobal

ggsave("outputs/figures/supp/SFig_sensi_phylogenetic_uncertainty_SMA.png",
       height = 5, width = 12.5, units = "in",
       plot = gglobal)

# PGLS-----
est_pgls <- data.frame(
  estimate = c(sensi$and_pgls$estimate, 
               sensi$gyn_pgls$estimate,
               sensi$pet_pgls$estimate,
               sensi$sep_pgls$estimate),
  organ = rep(c("and", "gyn", "pet", "sep"), each = nrow(sensi$and_physma)))


g3 <-
  ggplot(est_pgls %>% filter(organ %in% c("and", "gyn")), aes(x = estimate, fill = organ)) +
  geom_histogram(position = position_dodge(), 
                 bins = 20, show.legend = TRUE,
                 alpha = .8) +
  scale_fill_manual(values = cols[1:2], name = "", labels = c("Male", "Female")) +
  scale_x_continuous(breaks = seq(0.8, 1.3, 0.05)) +
  labs(y = "Frequency", x = "Estimated slopes",
       subtitle = "PGLS | (Male vs Female)") +
  tema(base_size = 20) +
  theme(legend.position = c(.88,.9),
        legend.background = element_blank());g3

g4 <-
  ggplot(est_pgls %>% filter(organ %in% c("pet", "sep")),
         aes(x = estimate, fill = organ)) +
  geom_histogram(position = position_dodge(), 
                 bins = 20, show.legend = TRUE,
                 alpha = .8) +
  scale_fill_manual(values = cols[3:4], name = "", labels = c("Petals", "Sepals")) +
  scale_x_continuous(breaks = seq(0.7, 1.5, 0.1)) +
  labs(y = "Frequency", x = "Estimated slopes",
       subtitle = "PGLS | (Petals vs Sepals)") +
  tema(base_size = 20) +
  theme(legend.position = c(.88,.9),
        legend.background = element_blank());g4

gpgls <- plot_grid(g3, g4,labels = LETTERS[1:2], 
                  label_size = 20, label_fontface = "plain"); gpgls

title <- ggdraw() + 
  draw_label("Sensitivity analysis - Phylogenetic uncertainty (PGLS)", size = 18)

gglobal <- 
  plot_grid(title, gpgls, ncol = 1, rel_heights = c(0.1, 1));gglobal

ggsave("outputs/figures/supp/SFig_sensi_phylogenetic_uncertainty_PGLS.png",
       height = 5, width = 12.5, units = "in",
       plot = gglobal)

### END------------