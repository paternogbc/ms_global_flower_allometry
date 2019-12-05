# Supp. Mat - Sensitivity analysis (axis independence)
# Allometric scalling with axis independence (part ~ (total - part))

### Packages:-------------------------------------------------------------------
library(sensiPhy)
library(phylolm)
library(rr2)
library(tidyverse)
library(smatr)
library(phytools)

### Start-----------------------------------------------------------------------
rm(list = ls())
source("R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
d <- read.csv(file = "data/processed/data_flower_biomass_partition.csv")
rownames(d) <- d$sp
t <- read.tree(file = "data/processed/working_study_phylo_tree.tre")

# Match data and phylogenetic tree----------------------------------------------
d  <- sensiPhy::match_dataphy(sp ~ 1, data = d, phy = t)

### Allometric regressions (Axis independence)-------------
### 1.1 OLS Male ----
olsand <- lm(log10(and) ~ log10(tot-and), data = d$data)
summary(olsand)
write.csv(summary(olsand)[[4]], file = "outputs/temp/model_estimates/sensi_axis_and_ols.csv")

### 1.2 PGLS Male ----
phyand <- phylolm(log10(and) ~ log10(tot - and), data = d$data, phy = d$phy,
                  model = "BM",
                  boot = 1000)
summary(phyand)
write.csv(summary(phyand)[[2]], file = "outputs/temp/model_estimates/sensi_axis_and_pgls.csv")

### 1.3 RMA Male ----
rmaand <- sma(log10(and) ~ log10(tot - and), data = d$data, slope.test = 1)
summary(rmaand)
write.csv(rmaand$groupsummary, file = "outputs/temp/model_estimates/sensi_axis_and_rma.csv")

### 1.4 RMA.phy Male ----
rmaphyand <- phyl.RMA(x = log10(extract(x = d$data, "tot") -
                                  extract(x = d$data, "and")), 
                      y = log10(extract(x = d$data, "and")), 
                      tree = d$phy, method = "BM")
summary.phyl.RMA(rmaphyand)
write.csv(summary.phyl.RMA(rmaphyand), file = "outputs/temp/model_estimates/sensi_axis_and_phyrma.csv")

### 2.1 OLS Female ----
olsgyn <- lm(log10(gyn) ~ log10(tot - gyn), data = d$data)
summary(olsgyn)
write.csv(summary(olsgyn)[[4]], file = "outputs/temp/model_estimates/sensi_axis_gyn_ols.csv")

### 2.2 PGLS Female ----
phygyn <- phylolm(log10(gyn) ~ log10(tot - gyn), data = d$data, phy = d$phy,
                  model = "BM",
                  boot = 1000)
summary(phygyn)
write.csv(summary(phygyn)[[2]], file = "outputs/temp/model_estimates/sensi_axis_gyn_pgls.csv")

### 2.3 RMA Female ----
rmagyn <- sma(log10(gyn) ~ log10(tot - gyn), data = d$data, slope.test = 1)
summary(rmagyn)
write.csv(rmagyn$groupsummary, file = "outputs/temp/model_estimates/sensi_axis_gyn_rma.csv")

### 2.4 RMA.phy Feamle ----
rmaphygyn <- phyl.RMA(x = log10(extract(x = d$data, "tot") -
                                  extract(x = d$data, "gyn")), 
                      y = log10(extract(x = d$data, "gyn")), 
                      tree = d$phy, method = "BM")
summary.phyl.RMA(rmaphygyn)
write.csv(summary.phyl.RMA(rmaphygyn), file = "outputs/temp/model_estimates/sensi_axis_gyn_phyrma.csv")

### 3.1 OLS Petals ----
olspet <- lm(log10(pet) ~ log10(tot - pet), data = d$data)
summary(olspet)
write.csv(summary(olspet)[[4]], file = "outputs/temp/model_estimates/sensi_axis_pet_ols.csv")

### 3.2 PGLS Petals ----
phypet <- phylolm(log10(pet) ~ log10(tot - pet), data = d$data, phy = d$phy,
                  model = "BM",
                  boot = 1000)
summary(phypet)
write.csv(summary(phypet)[[2]], file = "outputs/temp/model_estimates/sensi_axis_pet_pgls.csv")

### 3.3 RMA Petals ----
rmapet <- sma(log10(pet) ~ log10(tot - pet), data = d$data, slope.test = 1)
summary(rmapet)
write.csv(rmapet$groupsummary, file = "outputs/temp/model_estimates/sensi_axis_pet_rma.csv")

### 3.4 RMA.phy Petals ----
rmaphypet <- phyl.RMA(x = log10(extract(x = d$data, "tot") -
                                  extract(x = d$data, "pet")), 
                      y = log10(extract(x = d$data, "pet")), 
                      tree = d$phy, method = "BM")
summary.phyl.RMA(rmaphypet)
write.csv(summary.phyl.RMA(rmaphypet), file = "outputs/temp/model_estimates/sensi_axis_pet_phyrma.csv")

### 4.1 OLS Sepals ----
olssep <- lm(log10(sep) ~ log10(tot - sep), data = d$data)
summary(olssep)
write.csv(summary(olssep)[[4]], file = "outputs/temp/model_estimates/sensi_axis_sep_ols.csv")

### 4.2 PGLS Sepals ----
physep <- phylolm(log10(sep) ~ log10(tot - sep), data = d$data, phy = d$phy,
                  model = "BM",
                  boot = 1000)
summary(physep)
write.csv(summary(physep)[[2]], file = "outputs/temp/model_estimates/sensi_axis_sep_pgls.csv")

### 4.3 RMA Sepals ----
rmasep <- sma(log10(sep) ~ log10(tot - sep), data = d$data, slope.test = 1)
summary(rmasep)
write.csv(rmasep$groupsummary, file = "outputs/temp/model_estimates/sensi_axis_sep_rma.csv")

### 4.4 RMA.phy Sepals ----
rmaphysep <- phyl.RMA(x = log10(extract(x = d$data, "tot") - 
                                  extract(x = d$data, "sep")), 
                      y = log10(extract(x = d$data, "sep")), 
                      tree = d$phy, method = "BM")
summary.phyl.RMA(rmaphysep)
write.csv(summary.phyl.RMA(rmaphysep), file = "outputs/temp/model_estimates/sensi_axis_sep_phyrma.csv")

### 5. Prepare allometric scalling Tables----
### 5.1 OLS-----
table.ols <- 
  tab.ols(mod.list = list(olsand, olsgyn, olspet, olssep),
          title = "OLS regression (y ~ [x - Y])", 
          model.names = c("Male ~ Flower", 
                          "Female ~ Flower",
                          "Petals ~ Flower",
                          "Sepals ~ Flower"))
table.ols

### 5.2 PGLS-----
table.pgls <- 
  tab.phylolm(mod.list = list(phyand, phygyn, phypet, physep),
              phy = d$phy, title = "PGLS regression (y ~ [x - Y])", 
              model.names = c("Male ~ Flower", 
                              "Female ~ Flower",
                              "Petals ~ Flower",
                              "Sepals ~ Flower"))
table.pgls

### 5.3 RMA-----
table.rma <- 
  tab.rma(mod.list = list(rmaand, rmagyn, rmapet, rmasep),
          title = "RMA regression (y ~ [x - Y])", 
          model.names = c("Male ~ Flower", 
                          "Female ~ Flower",
                          "Petals ~ Flower",
                          "Sepals ~ Flower"))
table.rma

### 5.4 RMA.phy-----
table.phy.rma <- 
  tab.phy.rma(mod.list = list(rmaphyand, rmaphygyn, rmaphypet, rmaphysep),
              title = "RMA.phy regression (y ~ [x - Y])", 
              model.names = c("Male ~ Flower", 
                              "Female ~ Flower",
                              "Petals ~ Flower",
                              "Sepals ~ Flower"))
table.phy.rma

### joined table----
tab.scalling <- rbind(table.ols, table.pgls, table.rma, table.phy.rma)
tab.scalling

### * Save Model statistics----
write_csv(tab.scalling, "outputs/tables/supp/STable_sensi_axis_model_stats_allometric_scalling.csv")

### * Save fitted models----
models <- list(olsand,
               olsgyn,
               olspet,
               olssep,
               phyand,
               phygyn,
               phypet,
               physep,
               rmaand,
               rmagyn,
               rmapet,
               rmasep,
               rmaphyand,
               rmaphygyn,
               rmaphypet,
               rmaphysep)

names(models) <- c("olsand",
                   "olsgyn",
                   "olspet",
                   "olssep",
                   "phyand",
                   "phygyn",
                   "phypet",
                   "physep",
                   "rmaand",
                   "rmagyn",
                   "rmapet",
                   "rmasep",
                   "rmaphyand",
                   "rmaphygyn",
                   "rmaphypet",
                   "rmaphysep")
saveRDS(object = models, 
        file = "outputs/temp/sensi_axis_fitted_models_allometric_scaling.Rds")
### END-----  