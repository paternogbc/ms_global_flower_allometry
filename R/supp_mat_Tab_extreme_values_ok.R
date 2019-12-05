# Supp. Mat
# Table with extreme values in flower dimentions

### Packages:-------------------------------------------------------------------
library(tidyverse)

### Start-----------------------------------------------------------------------
rm(list = ls())
source("R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
d <- read.csv(file = "data/processed/data_flower_biomass_partition.csv")

### Calculate flower malenes
dc <- rename(d, 
            Flower = tot,
            Species = sp,
            Androecium = and,
            Gynoecium = gyn,
            Petals = pet,
            Sepals = sep) %>% 
  select(Species, Flower, everything()) %>% 
  mutate(maleness = Androecium / Gynoecium)
  
# Top 10 larger flowers---------------------------------------------------------
dc %>% 
  top_n(n = 10, wt = Flower) %>% 
  arrange(desc(Flower)) %>% 
  mutate(Rank = seq_len(10)) %>% 
  select(Rank, everything()) -> top
top

# Top 10 smaller flowers---------------------------------------------------------
dc %>% 
  top_n(n = 10, wt = desc(Flower)) %>% 
  arrange(Flower) %>% 
  mutate(Rank = 307:298) %>% 
  select(Rank, everything()) -> low
low

extremes <- rbind(top, low) 
extremes

# Save tables-------------------------------------------------------------------
write_csv(extremes, "outputs/tables/supp/STable_extreme_values_flower_biomass.csv")
