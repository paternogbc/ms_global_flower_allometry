# Supp. Mat
# test difference between allometric slopes
# Warning: Long time running (~1 hour)

# Packages----------------------------------------------------------------------
rm(list = ls())
library(sensiPhy)
library(phytools)
library(tidyverse)
library(smatr)
library(reshape2)
source("R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
d     <- read.csv("data/processed/data_flower_biomass_partition.csv")
p     <- read.tree("data/processed/working_study_phylo_tree.tre")

# Match data and tree
rownames(d) <- d$sp
dp <- match_dataphy(formula = sp ~ 1, data = d, phy = p)

# Reshape data
dmg <- dplyr::select(d, and, gyn, tot, sp)
dmg <- melt(dmg, id.vars = c("tot", "sp"), variable.name = "organ")
dps <- dplyr::select(d, pet, sep, tot, sp)
dps <- melt(dps, id.vars = c("tot", "sp"), variable.name = "organ")

### Test for difference in slopes between flower organs----
### 1 SMA---------------------
### 1.1 Male vs Female--------
m1   <- sma(log10(value) ~ log10(tot) * organ, data = dmg)
test <- m1$commoncoef
tab1.1 <- data.frame(method = "SMA", 
                   Slope_contrast = "Male vs Female",
                   LR = test$LR, df = test$df, P = test$p)
tab1.1
### 1.2 Petals vs Sepals----
m2   <- sma(log10(value) ~ log10(tot) * organ, data = dps)
test <- m2$commoncoef
tab1.2 <- data.frame(method = "SMA", 
                   Slope_contrast = "Petals vs Sepals",
                   LR = test$LR, df = test$df, P = test$p)

tab_sma <- rbind(tab1.1, tab1.2);tab_sma

### 2 phySMA---------------------
# (Permutation test)
set.seed(12345)

### 2.1 Male vs Female----------------------------------------------------------
test21 <- slope_test_phySMA(data = d, phy = p, groupvar = c("and","gyn"), xvar = "tot",
                           n_sim = 10000, method = "BM", type = "one-tailed")

### 2.2 Petals vs. Sepals-------------------------------------------------------
test22 <- slope_test_phySMA(data = d, phy = p, groupvar = c("pet","sep"), xvar = "tot",
                            n_sim = 10000, method = "BM", type = "one-tailed")

tab_physma <- rbind(test21$slope_test, test22$slope_test);tab_physma

### 3. Save tables--------------------------------------------------------------
write.csv(x = tab_sma, file = "outputs/tables/supp/STable_common_slope_test_SMA.csv")
write.csv(x = tab_physma, file = "outputs/tables/supp/STable_common_slope_test_phySMA.csv")

### Save simulation cache
saveRDS(object = list(test21, test22), file = "outputs/temp/common_slope_test_cache_phySMA.Rds")
### END----
