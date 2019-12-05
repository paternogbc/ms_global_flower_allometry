# Supp. Mat - Sensitivity analysis (Showy organs)

### Packages:-------------------------------------------------------------------
library(sensiPhy)
library(phylolm)
library(tidyverse)
library(smatr)
library(phytools)
library(cowplot)
library(rr2)

### Start-----------------------------------------------------------------------
rm(list = ls())
source("R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
p   <- read.tree("data/processed/working_study_phylo_tree.tre")
d   <- read.csv("data/processed/data_flower_biomass_partition.csv")
s   <- read.csv("data/raw/data_species_flower_organs_sensi.csv")
estimates <- read.csv("outputs/tables/supp/STable_model_stats_allometric_scalling.csv") 

cols <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")

# Merge data with taxonomic table
d <- left_join(d, s) %>%  as.data.frame()

d %>% 
  mutate(showy_stamens = ifelse(showy_stamens == "", yes = FALSE, no = TRUE),
         showy_sepals = ifelse(showy_sepals == "", yes = FALSE, no = TRUE),
         tepals = ifelse(tepals == "", yes = FALSE, no = TRUE),
         inflorecence = ifelse(inflorecence == "", yes = FALSE, no = TRUE)) -> d

# Match with Phylogeny----------------------------------------------------------
rownames(d) <- d$sp
d <- match_dataphy(tot ~ 1, data = d, phy = p)

### Showy organs, Tepals & inflorecense----
sp_sta <- d$data %>% filter(showy_stamens == FALSE) %>% pull(sp)
sp_sep <- d$data %>% filter(showy_sepals == FALSE) %>% pull(sp)
sp_tep <- d$data %>% filter(tepals == FALSE) %>% pull(sp)
sp_inf <- d$data %>% filter(inflorecence == FALSE) %>% pull(sp)

sp_removal <- list("showy stamens" = sp_sta,
                   "showy sepals" = sp_sep,
                   tepals = sp_tep,
                   inflorescence = sp_inf)

model_names <- c("Male ~ Flower", 
                 "Female ~ Flower",
                 "Petals ~ Flower",
                 "Sepals ~ Flower")

### SMA----
sma_tab <- data.frame()
for (i in seq_len(length(sp_removal))) {
  remain   <- sp_removal[[i]]
  name     <- names(sp_removal)[i]
  
  dc <- d$data %>% filter(sp %in% remain)
  
  and <- sma(log10(and) ~ log10(tot), data = dc, slope.test = 1)
  gyn <- sma(log10(gyn) ~ log10(tot), data = dc, slope.test = 1)
  pet <- sma(log10(pet) ~ log10(tot), data = dc, slope.test = 1)
  sep <- sma(log10(sep) ~ log10(tot), data = dc, slope.test = 1)
  
  tab <- tab.rma(mod.list = list(and, gyn, pet, sep), title = name, 
                 model.names = model_names)
  sma_tab <- rbind(sma_tab, tab)
}

sma_tab$Subtitle <- factor(sma_tab$Subtitle, levels = model_names)


### phySMA----
physma_tab <- data.frame()

for (i in seq_len(length(sp_removal))) {
  remain   <- sp_removal[[i]]
  name     <- names(sp_removal)[i]
  
  dc <- d$data %>% filter(sp %in% remain)
  rownames(dc) <- dc$sp
  dc <- match_dataphy(tot ~ 1, data = dc, phy = d$phy)
  and <- phyl.RMA(x = log10(extract(x = dc$data, "tot")), 
                  y = log10(extract(x = dc$data, "and")), 
                  tree = dc$phy, method = "BM")
  gyn <- phyl.RMA(x = log10(extract(x = dc$data, "tot")), 
                  y = log10(extract(x = dc$data, "gyn")), 
                  tree = dc$phy, method = "BM")
  pet <- phyl.RMA(x = log10(extract(x = dc$data, "tot")), 
                  y = log10(extract(x = dc$data, "pet")), 
                  tree = dc$phy, method = "BM")
  sep <- phyl.RMA(x = log10(extract(x = dc$data, "tot")), 
                  y = log10(extract(x = dc$data, "sep")), 
                  tree = dc$phy, method = "BM")
  
  tab <- tab.phy.rma(mod.list = list(and, gyn, pet, sep), title = name, 
                     model.names = model_names)
  physma_tab <- rbind(physma_tab, tab)
}
physma_tab$Subtitle <- factor(physma_tab$Subtitle, levels = model_names)

stat_tab1 <- rbind(sma_tab, physma_tab)
stat_tab1$model <- rep(c("SMA", "phySMA"), each = nrow(sma_tab))


### OLS----
ols_tab <- data.frame()
for (i in seq_len(length(sp_removal))) {
  remain   <- sp_removal[[i]]
  name     <- names(sp_removal)[i]
  
  dc <- d$data %>% filter(sp %in% remain)
  
  and <- sma(log10(and) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  gyn <- sma(log10(gyn) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  pet <- sma(log10(pet) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  sep <- sma(log10(sep) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  
  tab <- tab.rma(mod.list = list(and, gyn, pet, sep), title = name, 
                 model.names = model_names)
  ols_tab <- rbind(ols_tab, tab)
}

ols_tab$Subtitle <- factor(ols_tab$Subtitle, levels = model_names)

### PGLS----
pgls_tab <- data.frame()
for (i in seq_len(length(sp_removal))) {
  remain   <- sp_removal[[i]]
  name     <- names(sp_removal)[i]
  
  dc <- d$data %>% filter(sp %in% remain)
  rownames(dc) <- dc$sp
  dc <- match_dataphy(tot ~ 1, data = dc, phy = d$phy)
  
  and <- phylolm(log10(and) ~ log10(tot), data = dc$data, phy = dc$phy,
                 model = "BM",
                 boot = 1000)
  gyn <- phylolm(log10(gyn) ~ log10(tot), data = dc$data, phy = dc$phy,
                 model = "BM",
                 boot = 1000)
  pet <- phylolm(log10(pet) ~ log10(tot), data = dc$data, phy = dc$phy,
                 model = "BM",
                 boot = 1000)
  sep <- phylolm(log10(sep) ~ log10(tot), data = dc$data, phy = dc$phy,
                 model = "BM",
                 boot = 1000)
  
  tab <- tab.phylolm(mod.list = list(and, gyn, pet, sep), title = name, 
                 model.names = model_names, phy = dc$phy)
  pgls_tab <- rbind(pgls_tab, tab)
}

pgls_tab$Subtitle <- factor(pgls_tab$Subtitle, levels = model_names)

stat_tab2 <- rbind(ols_tab, pgls_tab)
stat_tab2$model <- rep(c("OLS", "PGLS"), each = nrow(ols_tab))

# Figures---------------
# SMA & phySMA-----
es <- estimates %>% filter(Title %in% c("RMA regression", "RMA.phy regression"))
vline.dat <- data.frame(Subtitle = model_names, 
                        slo = es$slo,
                        model = rep(c("SMA", "phySMA"), each = 4)) 

stat_tab1$model <- factor(stat_tab1$model, levels = c( "SMA", "phySMA"))

g1 <- ggplot(stat_tab1, aes(y = slo, x = Title, 
                           color = Subtitle, 
                           ymin = slo.lo, ymax = slo.hi)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_hline(data = vline.dat, aes(yintercept = slo, color = Subtitle),
             show.legend = F, size = 2, alpha = .8) +
  geom_point(show.legend = F, size = 6) +
  geom_errorbar(show.legend = FALSE, width = .2) +
  scale_color_manual(values = cols[c(2,1,3,4)]) +
  facet_grid(model ~ Subtitle) +
  labs(x = "Species set removed", y = "Estimated slope") +
  tema(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90));g1

title <- ggdraw() + draw_label("Sensitivity analysis - Species variation in showiness strategy", size = 20)
  
gglobal <- plot_grid(title, g1, ncol = 1, rel_heights = c(0.1, 1));gglobal
  
### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/supp/SFig_sensi_species_set_SMA.png",
         height = 8.5, width = 11.2, units = "in",
         plot = gglobal)
  
# OLS & PGLS-----
es <- estimates %>% filter(Title %in% c("OLS regression", "PGLS regression"))
vline.dat <- data.frame(Subtitle = model_names, 
                        slo = es$slo,
                        model = rep(c("OLS", "PGLS"), each = 4)) 

stat_tab2$model <- factor(stat_tab2$model, levels = c( "OLS", "PGLS"))

g1 <- ggplot(stat_tab2, aes(y = slo, x = Title, 
                            color = Subtitle, 
                            ymin = slo.lo, ymax = slo.hi)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_hline(data = vline.dat, aes(yintercept = slo, color = Subtitle),
             show.legend = F, size = 2, alpha = .8) +
  geom_point(show.legend = F, size = 6) +
  geom_errorbar(show.legend = FALSE, width = .2) +
  scale_color_manual(values = cols[c(2,1,3,4)]) +
  facet_grid(model ~ Subtitle) +
  labs(x = "Species set removed", y = "Estimated slope") +
  tema(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90));g1

title <- ggdraw() + draw_label("Sensitivity analysis - Species variation in showiness strategy", size = 20)

gglobal <- plot_grid(title, g1, ncol = 1, rel_heights = c(0.1, 1));gglobal

### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/supp/SFig_sensi_species_set_OLS.png",
       height = 8.5, width = 11.2, units = "in",
       plot = gglobal)

# P-value visualization (just for checking B == 1)
stat_tab <- rbind(stat_tab1, stat_tab2)  %>% 
  select(species_removed = Title, fit = Subtitle, model, everything())

stat_tab %>% 
  ggplot(aes(y = p_b1, x = species_removed, fill = fit)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0.05) +
  scale_fill_manual(values = cols) +
  facet_grid(model ~.) +
  tema(base_size = 18)

### Save tables---------------
write_csv(stat_tab, 
          path = "outputs/tables/supp/STable_sensi_species_set.csv")
### END----