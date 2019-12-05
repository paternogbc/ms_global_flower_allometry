# Supp. Mat
# Build higher level taxonomic table

### Packages:-------------------------------------------------------------------
library(tidyverse)
library(Taxonstand)
library(stdnames) # Available at: https://github.com/paternogbc/gazp-toolbox/tree/master/code/stdnames

# Load functions ----------------------------------------------------------
source(file = "R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
d  <- read_csv("data/processed/data_flower_biomass_partition.csv")

# 1. Standardize species names after The Plant List-----------------------------
#out <- std_names(x = data.frame(d), species_column = "sp", version = "1.1")
#saveRDS(object = out, file = "outputs/temp/output_plant_names_standardization.Rds")
out <- readRDS(file = "outputs/temp/output_plant_names_standardization.Rds")

# Save standardized plant names
taxonomy <- out$detailed_output$tpl_full
taxonomy <-
  taxonomy %>% 
  select(order, family, original_name, tpl_genus, tpl_name, 
         tpl_authority, tpl_id) %>% 
  arrange(order, family, tpl_genus, tpl_name)

# Fix family to match the plant list: Richea_continenti -> Ericaceae
taxonomy[taxonomy$original_name == "Richea_continentis", ]$family = "Ericaceae"
taxonomy[taxonomy$original_name == "Richea_continentis", ]$order = "Ericales"

head(taxonomy)

# Save table with higher taxonomy
write_csv(taxonomy, "outputs/tables/supp/STable_plant_names_taxonomy.csv")
