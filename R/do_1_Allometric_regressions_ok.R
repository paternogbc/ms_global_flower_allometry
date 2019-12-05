### Allometric scalling of flower biomass partition

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

### Allometric regressions-------------
### 1.1 OLS Male ----
olsand <- lm(log10(and) ~ log10(tot), data = d$data)
summary(olsand)
write.csv(summary(olsand)[[4]], file = "outputs/temp/model_estimates/and_ols.csv")

### 1.2 PGLS Male ----
phyand <- phylolm(log10(and) ~ log10(tot), data = d$data, phy = d$phy,
                  model = "BM",
                  boot = 1000)
summary(phyand)
write.csv(summary(phyand)[[2]], file = "outputs/temp/model_estimates/and_pgls.csv")

### Fit with lambda
phyand2 <- phylolm(log10(and) ~ log10(tot), data = d$data, phy = d$phy,
                  model = "lambda", boot = 100)
summary(phyand2)

### 1.3 RMA Male ----
rmaand <- sma(log10(and) ~ log10(tot), data = d$data, slope.test = 1)
summary(rmaand)
write.csv(rmaand$groupsummary, file = "outputs/temp/model_estimates/and_rma.csv")

### 1.4 RMA.phy Male ----
rmaphyand <- phyl.RMA(x = log10(extract(x = d$data, "tot")), 
                      y = log10(extract(x = d$data, "and")), 
                      tree = d$phy, method = "BM")
summary.phyl.RMA(rmaphyand)
write.csv(summary.phyl.RMA(rmaphyand), file = "outputs/temp/model_estimates/and_phyrma.csv")

### Fit with lambda
rmaphyand2 <- phyl.RMA(x = log10(extract(x = d$data, "tot")), 
                      y = log10(extract(x = d$data, "and")), 
                      tree = d$phy, method = "lambda")
summary.phyl.RMA(rmaphyand2)

### 2.1 OLS Female ----
olsgyn <- lm(log10(gyn) ~ log10(tot), data = d$data)
summary(olsgyn)
write.csv(summary(olsgyn)[[4]], file = "outputs/temp/model_estimates/gyn_ols.csv")

### 2.2 PGLS Female ----
phygyn <- phylolm(log10(gyn) ~ log10(tot), data = d$data, phy = d$phy,
                  model = "BM",
                  boot = 1000)
summary(phygyn)
write.csv(summary(phygyn)[[2]], file = "outputs/temp/model_estimates/gyn_pgls.csv")

### Fit with lambda
phygyn2 <- phylolm(log10(gyn) ~ log10(tot), data = d$data, phy = d$phy,
                   model = "lambda", boot = 100)
summary(phygyn2)

### 2.3 RMA Female ----
rmagyn <- sma(log10(gyn) ~ log10(tot), data = d$data, slope.test = 1)
summary(rmagyn)
write.csv(rmagyn$groupsummary, file = "outputs/temp/model_estimates/gyn_rma.csv")

### 2.4 RMA.phy Feamle ----
rmaphygyn <- phyl.RMA(x = log10(extract(x = d$data, "tot")), 
                      y = log10(extract(x = d$data, "gyn")), 
                      tree = d$phy, method = "BM")
summary.phyl.RMA(rmaphygyn)
write.csv(summary.phyl.RMA(rmaphygyn), file = "outputs/temp/model_estimates/gyn_phyrma.csv")

# Fit with lambda
rmaphygyn2 <- phyl.RMA(x = log10(extract(x = d$data, "tot")), 
                      y = log10(extract(x = d$data, "gyn")), 
                      tree = d$phy, method = "lambda")
summary.phyl.RMA(rmaphygyn2)

### 3.1 OLS Petals ----
olspet <- lm(log10(pet) ~ log10(tot), data = d$data)
summary(olspet)
write.csv(summary(olspet)[[4]], file = "outputs/temp/model_estimates/pet_ols.csv")

### 3.2 PGLS Petals ----
phypet <- phylolm(log10(pet) ~ log10(tot), data = d$data, phy = d$phy,
                  model = "BM",
                  boot = 1000)
summary(phypet)
write.csv(summary(phypet)[[2]], file = "outputs/temp/model_estimates/pet_pgls.csv")

# Fitted with lambda
phypet2 <- phylolm(log10(pet) ~ log10(tot), data = d$data, phy = d$phy,
                  model = "lambda", boot = 100)
summary(phypet2)

### 3.3 RMA Petals ----
rmapet <- sma(log10(pet) ~ log10(tot), data = d$data, slope.test = 1)
summary(rmapet)
write.csv(rmapet$groupsummary, file = "outputs/temp/model_estimates/pet_rma.csv")

### 3.4 RMA.phy Petals ----
rmaphypet <- phyl.RMA(x = log10(extract(x = d$data, "tot")), 
                      y = log10(extract(x = d$data, "pet")), 
                      tree = d$phy, method = "BM")
summary.phyl.RMA(rmaphypet)
write.csv(summary.phyl.RMA(rmaphypet), file = "outputs/temp/model_estimates/pet_phyrma.csv")

# Fitted with lambda
rmaphypet2 <- phyl.RMA(x = log10(extract(x = d$data, "tot")), 
                      y = log10(extract(x = d$data, "pet")), 
                      tree = d$phy, method = "lambda")
summary.phyl.RMA(rmaphypet2)

### 4.1 OLS Sepals ----
olssep <- lm(log10(sep) ~ log10(tot), data = d$data)
summary(olssep)
write.csv(summary(olssep)[[4]], file = "outputs/temp/model_estimates/sep_ols.csv")

### 4.2 PGLS Sepals ----
physep <- phylolm(log10(sep) ~ log10(tot), data = d$data, phy = d$phy,
                  model = "BM",
                  boot = 1000)
summary(physep)
write.csv(summary(physep)[[2]], file = "outputs/temp/model_estimates/sep_pgls.csv")

# Fitted with lambda
physep2 <- phylolm(log10(sep) ~ log10(tot), data = d$data, phy = d$phy,
                  model = "lambda", boot = 100)
summary(physep2)

### 4.3 RMA Sepals ----
rmasep <- sma(log10(sep) ~ log10(tot), data = d$data, slope.test = 1)
summary(rmasep)
write.csv(rmasep$groupsummary, file = "outputs/temp/model_estimates/sep_rma.csv")

### 4.4 RMA.phy Sepals ----
rmaphysep <- phyl.RMA(x = log10(extract(x = d$data, "tot")), 
                      y = log10(extract(x = d$data, "sep")), 
                      tree = d$phy, method = "BM")
summary.phyl.RMA(rmaphysep)
write.csv(summary.phyl.RMA(rmaphysep), file = "outputs/temp/model_estimates/sep_phyrma.csv")

# Fitted with lambda
rmaphysep2 <- phyl.RMA(x = log10(extract(x = d$data, "tot")), 
                      y = log10(extract(x = d$data, "sep")), 
                      tree = d$phy, method = "lambda")
summary.phyl.RMA(rmaphysep2)

### 5. Prepare allometric scalling Tables----
### 5.1 OLS-----
table.ols <- 
  tab.ols(mod.list = list(olsand, olsgyn, olspet, olssep),
          title = "OLS regression", 
          model.names = c("Male ~ Flower", 
                          "Female ~ Flower",
                          "Petals ~ Flower",
                          "Sepals ~ Flower"))
table.ols

### 5.2 PGLS-----
table.pgls <- 
  tab.phylolm(mod.list = list(phyand, phygyn, phypet, physep),
              phy = d$phy, title = "PGLS regression", 
              model.names = c("Male ~ Flower", 
                              "Female ~ Flower",
                              "Petals ~ Flower",
                              "Sepals ~ Flower"))
table.pgls

# Fitted with lambda
table_pgls_lam <-
  tab.phylolm(mod.list = list(phyand2, phygyn2, phypet2, physep2),
              phy = d$phy, title = "PGLS regression - lambda", 
              model.names = c("Male ~ Flower", 
                              "Female ~ Flower",
                              "Petals ~ Flower",
                              "Sepals ~ Flower"))
table_pgls_lam$lambda <- 
  c(phyand2$optpar, phygyn2$optpar, phypet2$optpar, physep2$optpar)

### 5.3 RMA-----
table.rma <- 
  tab.rma(mod.list = list(rmaand, rmagyn, rmapet, rmasep),
          title = "RMA regression", 
          model.names = c("Male ~ Flower", 
                          "Female ~ Flower",
                          "Petals ~ Flower",
                          "Sepals ~ Flower"))
table.rma

### 5.4 RMA.phy-----
table.phy.rma <- 
  tab.phy.rma(mod.list = list(rmaphyand, rmaphygyn, rmaphypet, rmaphysep),
              title = "RMA.phy regression", 
              model.names = c("Male ~ Flower", 
                              "Female ~ Flower",
                              "Petals ~ Flower",
                              "Sepals ~ Flower"))
table.phy.rma

# Fitted with lambda
table.phy.rma_lam <- 
  tab.phy.rma(mod.list = list(rmaphyand2, rmaphygyn2, rmaphypet2, rmaphysep2),
              title = "RMA.phy regression", 
              model.names = c("Male ~ Flower", 
                              "Female ~ Flower",
                              "Petals ~ Flower",
                              "Sepals ~ Flower"))

table.phy.rma_lam$lambda <- 
  c(rmaphyand2$lambda, rmaphygyn2$lambda, rmaphypet2$lambda, rmaphysep2$lambda)

### joined table----
tab.scalling <- rbind(table.ols, table.pgls, table.rma, table.phy.rma)
tab.scalling

tab.scalling_lam <- rbind(table.ols, table.pgls, table.rma, table.phy.rma)
tab.scalling_lam <- rbind(table_pgls_lam, table.phy.rma_lam)

### * Save Model statistics----
write_csv(tab.scalling, "outputs/tables/supp/STable_model_stats_allometric_scalling.csv")
write_csv(tab.scalling_lam, "outputs/tables/supp/STable_model_stats_allometric_scalling_lambda.csv")

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
        file = "outputs/temp/fitted_models_allometric_scaling.Rds")
### END-----  