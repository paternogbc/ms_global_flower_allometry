# Supplementary material
# Figure:Phylogenetic tree with orders colored

# Packages ----------------------------------------------------------------
library(ggtree)
library(RColorBrewer)
library(tidyverse)
library(phytools)
library(sensiPhy)

### Load data:-----
rm(list = ls())
source("R/ZZZ_functions.R")

dc <- read.csv("outputs/tables/supp/STable_plant_names_taxonomy.csv")
pp <- read.tree("data/processed/working_study_phylo_tree.tre")

rownames(dc) <- dc$original_name
d <- match_dataphy(formula = original_name ~ 1,
                   data = dc,
                   phy = pp)

### find species by Order
sp.ord <-
  split(as.character(d$data$original_name),
        as.character(d$data$order))
ord.na <- names(sp.ord)

top_ord <-
  dc %>%
  group_by(order) %>%
  count() %>%
  pull(order) %>%
  as.character()

### Find RMCA for each order:
ord.na <- ord.na[ord.na %in% top_ord]
nods <- list()

for (i in ord.na) {
  nod <- findMRCA(tree = pp, sp.ord[[i]])
  if (is.null(nod)) {
    nod <- NA
  }
  nods[i] <- nod
}

# Exclude NA nods
nods <- nods[!is.na(nods)]
length(nods)

### Color vector:
col_vector = c(
  "#8dd3c7",
  "#bebada",
  "#fb8072",
  "#80b1d3",
  "#fdb462",
  "#b3de69",
  "#ec7014",
  "#bc80bd",
  "#00441b",
  "#1f78b4",
  "#33a02c",
  "#fb9a99",
  "#e31a1c",
  "#ff7f00",
  "#f781bf",
  "#08306b",
  "#b2182b",
  "#2166ac",
  "#1b7837",
  "#de77ae",
  "#8c510a",
  "#01665e"
)
set.seed(1311)
col_vector <- col_vector[sample(1:22)]

tt <- ggtree(pp, layout = "circular", ladderize = FALSE)
tt
for (i in 1:length(nods)) {
  tt <-
    tt +
    geom_hilight(node = nods[[i]],
                 fill = col_vector[i],
                 alpha = .7) +
    geom_cladelabel(
      color = col_vector[i],
      node = nods[[i]],
      label = names(nods)[i],
      angle = "auto",
      offset = 1,
      alpha = 1,
      fontsize = 4,
      barsize = 1
    )
  tt
}
tt <- tt + xlim(0, 200)
tt
ggsave(
  filename = "outputs/figures/supp/SFig_phylogeny_orders_colored.png",
  width = 12,
  height = 8,
  plot = tt
)
# END-----