# Supp. Mat - Sensitivity analysis (Taxonomic influence)

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
tax <- read_csv("outputs/tables/supp/STable_plant_names_taxonomy.csv")
estimates <- read.csv("outputs/tables/supp/STable_model_stats_allometric_scalling.csv") 

cols <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")

# Merge data with taxonomic table
d <- left_join(d, select(tax, order, family, sp = original_name)) %>%  as.data.frame()

# Match with Phylogeny----------------------------------------------------------
rownames(d) <- d$sp
d <- match_dataphy(tot ~ 1, data = d, phy = p)

### Sensitivity to Taxonomic influence----
### 1. Remove top 5 Order--------------------
top.ord <- as.character(count(d$data, order, sort = T)$order[1:5])
model_names <- c("Male ~ Flower", 
                 "Female ~ Flower",
                 "Petals ~ Flower",
                 "Sepals ~ Flower")
### SMA----
sma_ord <- data.frame()
for (i in top.ord) {
  dc <- filter(d$data, order != i)
  and <- sma(log10(and) ~ log10(tot), data = dc, slope.test = 1)
  gyn <- sma(log10(gyn) ~ log10(tot), data = dc, slope.test = 1)
  pet <- sma(log10(pet) ~ log10(tot), data = dc, slope.test = 1)
  sep <- sma(log10(sep) ~ log10(tot), data = dc, slope.test = 1)
  
  tab <- tab.rma(mod.list = list(and, gyn, pet, sep), title = i, 
          model.names = model_names)
  sma_ord <- rbind(sma_ord, tab)
}
sma_ord$Subtitle <- factor(sma_ord$Subtitle, levels = model_names)

### phySMA----
physma_ord <- data.frame()
for (i in top.ord) {
  dc <- filter(d$data, order != i)
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
  
  tab <- tab.phy.rma(mod.list = list(and, gyn, pet, sep), title = i, 
                 model.names = model_names)
  physma_ord <- rbind(physma_ord, tab)
}
physma_ord$Subtitle <- factor(physma_ord$Subtitle, levels = model_names)

### OLS----
ols_ord <- data.frame()
for (i in top.ord) {
  dc <- filter(d$data, order != i)
  and <- sma(log10(and) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  gyn <- sma(log10(gyn) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  pet <- sma(log10(pet) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  sep <- sma(log10(sep) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  
  tab <- tab.rma(mod.list = list(and, gyn, pet, sep), title = i, 
                 model.names = model_names)
  ols_ord <- rbind(ols_ord, tab)
}
ols_ord$Subtitle <- factor(ols_ord$Subtitle, levels = model_names)

### PGLS----
pgls_ord <- data.frame()
for (i in top.ord) {
  dc <- filter(d$data, order != i)
  rownames(dc) <- dc$sp
  dc <- match_dataphy(tot ~ 1, data = dc, phy = d$phy)
  
  and <- phylolm(log10(and) ~ log10(tot), data = dc$data, phy = dc$phy, 
                 model = "BM", boot = 1000)
  gyn <- phylolm(log10(gyn) ~ log10(tot), data = dc$data, phy = dc$phy, 
                 model = "BM", boot = 1000)
  pet <- phylolm(log10(pet) ~ log10(tot), data = dc$data, phy = dc$phy, 
                 model = "BM", boot = 1000)
  sep <- phylolm(log10(sep) ~ log10(tot), data = dc$data, phy = dc$phy, 
                 model = "BM", boot = 1000)
  
  tab <- tab.phylolm(mod.list = list(and, gyn, pet, sep), title = i, phy = dc$phy,
                 model.names = model_names)
  pgls_ord <- rbind(pgls_ord, tab)
}

pgls_ord$Subtitle <- factor(pgls_ord$Subtitle, levels = model_names)

# Figures---------------
# (SMA & phySMA)------
stat_tab_ord <- rbind(sma_ord, physma_ord)
stat_tab_ord$model <- rep(c("SMA", "phySMA"), each = 20)
stat_tab_ord$model <- factor(stat_tab_ord$model, levels = c("SMA", "phySMA"))

es <- estimates %>% filter(Title %in% c("RMA regression", "RMA.phy regression"))
vline.dat <- data.frame(Subtitle = model_names, 
                        slo = es$slo,
                        model = rep(c("SMA", "phySMA"), each = 4)) 
gord <- 
  ggplot(stat_tab_ord, aes(y = slo, x = reorder(Title, n), 
                    color = Subtitle, 
                    ymin = slo.lo, ymax = slo.hi)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_hline(data = vline.dat, aes(yintercept = slo, color = Subtitle),
             show.legend = F, size = 2, alpha = .8) +
  geom_point(show.legend = F, size = 5) +
  geom_errorbar(show.legend = FALSE, width = .2) +
  scale_color_manual(values = cols[c(2,1,3,4)]) +
  facet_grid(model ~ Subtitle) +
  labs(x = "Removed order", y = "Estimated slope") +
  tema(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5));gord

title <- ggdraw() + 
  draw_label("Sensitivity analysis - Taxonomic influence (SMA & phySMA)", size = 20)

gglobal <- 
  plot_grid(title, gord, ncol=1, rel_heights=c(0.1, 1));gglobal

### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/supp/SFig_sensi_taxonomic_influ_order_SMA.png",
       height = 8.5, width = 11.2, units = "in",
       plot = gglobal)


# (OLS & PGLS)-----
stat_tab_ord <- rbind(ols_ord, pgls_ord)
stat_tab_ord$model <- rep(c("OLS", "PGLS"), each = 20)
stat_tab_ord$model <- factor(stat_tab_ord$model, levels = c("OLS", "PGLS"))

es <- estimates %>% filter(Title %in% c("OLS regression", "PGLS regression"))
vline.dat <- data.frame(Subtitle = model_names, 
                        slo = es$slo,
                        model = rep( c("OLS", "PGLS"), each = 4)) 
gord <- 
  ggplot(stat_tab_ord, aes(y = slo, x = reorder(Title, n), 
                           color = Subtitle, 
                           ymin = slo.lo, ymax = slo.hi)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_hline(data = vline.dat, aes(yintercept = slo, color = Subtitle),
             show.legend = F, size = 2, alpha = .8) +
  geom_point(show.legend = F, size = 5) +
  geom_errorbar(show.legend = FALSE, width = .2) +
  scale_color_manual(values = cols[c(2,1,3,4)]) +
  facet_grid(model ~ Subtitle) +
  labs(x = "Removed order", y = "Estimated slope") +
  tema(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5));gord

title <- ggdraw() + 
  draw_label("Sensitivity analysis - Taxonomic influence (OLS & PGLS)", size = 20)

gglobal <- 
  plot_grid(title, gord, ncol=1, rel_heights=c(0.1, 1));gglobal

### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/supp/SFig_sensi_taxonomic_influ_order_OLS.png",
       height = 8.5, width = 11.2, units = "in",
       plot = gglobal)

### 2. Remove top 5 Families--------------------
top.fam <- as.character(count(d$data, family, sort = T)$family[1:5])

### SMA----
sma_fam <- data.frame()
for (i in top.fam) {
  dc <- filter(d$data, family != i)
  and <- sma(log10(and) ~ log10(tot), data = dc, slope.test = 1)
  gyn <- sma(log10(gyn) ~ log10(tot), data = dc, slope.test = 1)
  pet <- sma(log10(pet) ~ log10(tot), data = dc, slope.test = 1)
  sep <- sma(log10(sep) ~ log10(tot), data = dc, slope.test = 1)
  
  tab <- tab.rma(mod.list = list(and, gyn, pet, sep), title = i, 
                 model.names = model_names)
  sma_fam <- rbind(sma_fam, tab)
}
sma_fam$Subtitle <- factor(sma_fam$Subtitle, levels = model_names)

### phySMA----
physma_fam <- data.frame()
for (i in top.fam) {
  dc <- filter(d$data, family != i)
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
  
  tab <- tab.phy.rma(mod.list = list(and, gyn, pet, sep), title = i, 
                     model.names = model_names)
  physma_fam <- rbind(physma_fam, tab)
}
physma_fam$Subtitle <- factor(physma_fam$Subtitle, levels = model_names)

### OLS----
ols_fam <- data.frame()
for (i in top.fam) {
  dc <- filter(d$data, family != i)
  and <- sma(log10(and) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  gyn <- sma(log10(gyn) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  pet <- sma(log10(pet) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  sep <- sma(log10(sep) ~ log10(tot), data = dc, slope.test = 1, method = "OLS")
  
  tab <- tab.rma(mod.list = list(and, gyn, pet, sep), title = i, 
                 model.names = model_names)
  ols_fam <- rbind(ols_fam, tab)
}
ols_fam$Subtitle <- factor(ols_fam$Subtitle, levels = model_names)

### PGLS----
pgls_fam <- data.frame()
for (i in top.fam) {
  dc <- filter(d$data, family != i)
  rownames(dc) <- dc$sp
  dc <- match_dataphy(tot ~ 1, data = dc, phy = d$phy)
  
  and <- phylolm(log10(and) ~ log10(tot), data = dc$data, phy = dc$phy, 
                 model = "BM", boot = 1000)
  gyn <- phylolm(log10(gyn) ~ log10(tot), data = dc$data, phy = dc$phy, 
                 model = "BM", boot = 1000)
  pet <- phylolm(log10(pet) ~ log10(tot), data = dc$data, phy = dc$phy, 
                 model = "BM", boot = 1000)
  sep <- phylolm(log10(sep) ~ log10(tot), data = dc$data, phy = dc$phy, 
                 model = "BM", boot = 1000)
  
  tab <- tab.phylolm(mod.list = list(and, gyn, pet, sep), title = i, phy = dc$phy,
                     model.names = model_names)
  pgls_fam <- rbind(pgls_fam, tab)
}

pgls_fam$Subtitle <- factor(pgls_fam$Subtitle, levels = model_names)

# Figures---------------
# (SMA & phySMA)------
stat_tab <- rbind(sma_fam, physma_fam)
stat_tab$model <- rep(c("SMA", "phySMA"), each = 20)
stat_tab$model <- factor(stat_tab$model, levels = c("SMA", "phySMA"))

es <- estimates %>% filter(Title %in% c("RMA regression", "RMA.phy regression"))
vline.dat <- data.frame(Subtitle = model_names, 
                        slo = es$slo,
                        model = rep(c("SMA", "phySMA"), each = 4)) 
gfam <- 
  ggplot(stat_tab, aes(y = slo, x = reorder(Title, n), 
                           color = Subtitle, 
                           ymin = slo.lo, ymax = slo.hi)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_hline(data = vline.dat, aes(yintercept = slo, color = Subtitle),
             show.legend = F, size = 2, alpha = .8) +
  geom_point(show.legend = F, size = 5) +
  geom_errorbar(show.legend = FALSE, width = .2) +
  scale_color_manual(values = cols[c(2,1,3,4)]) +
  facet_grid(model ~ Subtitle) +
  labs(x = "Removed family", y = "Estimated slope") +
  tema(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5));gfam

title <- ggdraw() + 
  draw_label("Sensitivity analysis - Taxonomic influence (SMA & phySMA)", size = 20)

gglobal <- 
  plot_grid(title, gfam, ncol=1, rel_heights=c(0.1, 1));gglobal

### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/supp/SFig_sensi_taxonomic_influ_family_SMA.png",
       height = 8.5, width = 11.2, units = "in",
       plot = gglobal)

# (OLS & PGLS)-----
stat_tab <- rbind(ols_fam, pgls_fam)
stat_tab$model <- rep(c("OLS", "PGLS"), each = 20)
stat_tab$model <- factor(stat_tab$model, levels = c("OLS", "PGLS"))

es <- estimates %>% filter(Title %in% c("OLS regression", "PGLS regression"))
vline.dat <- data.frame(Subtitle = model_names, 
                        slo = es$slo,
                        model = rep( c("OLS", "PGLS"), each = 4)) 
gfam <- 
  ggplot(stat_tab, aes(y = slo, x = reorder(Title, n), 
                           color = Subtitle, 
                           ymin = slo.lo, ymax = slo.hi)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_hline(data = vline.dat, aes(yintercept = slo, color = Subtitle),
             show.legend = F, size = 2, alpha = .8) +
  geom_point(show.legend = F, size = 5) +
  geom_errorbar(show.legend = FALSE, width = .2) +
  scale_color_manual(values = cols[c(2,1,3,4)]) +
  facet_grid(model ~ Subtitle) +
  labs(x = "Removed family", y = "Estimated slope") +
  tema(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5));gfam

title <- ggdraw() + 
  draw_label("Sensitivity analysis - Taxonomic influence (OLS & PGLS)", size = 20)

gglobal <- 
  plot_grid(title, gfam, ncol = 1, rel_heights = c(0.1, 1));gglobal

### Export graphs to A4 landscape format (pdf):
ggsave("outputs/figures/supp/SFig_sensi_taxonomic_influ_family_OLS.png",
       height = 8.5, width = 11.2, units = "in",
       plot = gglobal)

### Save tables---------------
stat_tab <- rbind(sma_ord, )
write_csv(stat_tab, 
          path = "outputs/tables/supp/STable_sensi_taxonomic_influence.csv")
