# Supplementary material
# Figure: Phylogentic tree with flower biomass in bars

# Packages ----------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(sensiPhy)
library(phytools)
source(file = "R/ZZZ_functions.R")

# Load data ---------------------------------------------------------------
# Permutation test for phySMA
d <- read.csv("data/processed/data_flower_biomass_partition.csv")
p <- read.tree("data/processed/working_study_phylo_tree.tre")

# Plot tree-----
rownames(d) <- d$sp
flower <- sqrt(extract(x = d, column = "tot"))

# Square root of the flower biomass (scale)
pdf(file = "outputs/figures/supp/SFig_phylo_tree_flower_biomass.pdf")
plotTree.wBars(tree = p, x = flower, type = "fan", color = gray(.6), 
               border = "red", width = .5)
dev.off()
# END----