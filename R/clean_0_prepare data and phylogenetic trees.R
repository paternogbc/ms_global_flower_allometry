### Script to add missing species to the phylogeny.

### Packages:-------------------------------------------------------------------
library(sensiPhy)
library(tidyverse)

### Start-----------------------------------------------------------------------
rm(list = ls())
set.seed(12345)

### Load data-------------------------------------------------------------------
d  <- read.csv("data/raw/data_flower_biomass_partition.csv")
p  <- read.tree("data/raw/ALLOTB.tre")

### 1. Match data and phy-------------------------------------------------------
rownames(d) <- d$sp
dat <- sensiPhy::match_dataphy(formula = sp ~ 1, data = d, phy = p)

### Find species missing in phylogeny-------------------------------------------
missing <- dat$dropped[dat$dropped %in% d$sp]
write.csv(cbind(missing), file = "outputs/temp/Table_sp_initial_missing_at_phylogeny.csv")

### Get genus from data and phylogeny
mg <- sapply(strsplit(missing, "_"), function(x) x[[1]][1])
mp <- sapply(strsplit(p$tip.label, "_"), function(x) x[[1]][1])

### Find genus not present in the phylogeny
miss.genus <- mg %in% mp
missing <- missing[miss.genus]
mg <- mg[miss.genus]

### Gather potential species for replacement (the same genus)-------------------
replace <- list()
for (i in 1:length(mg)) {
  wp <- which(mp == mg[i])
  genus <- p$tip.label[wp]
  replace[[i]] <- genus
}
names(replace) <- mg

### Count the number of potential species by genus
replace_n <- data.frame("genus" = character(length = length(replace)),
                        "N" = numeric(length = length(replace)),
                        stringsAsFactors = F)
for (j in 1:length(replace)) {
  g <- names(replace)[j]
  l <- length(replace[[j]])
  replace_n[j, ] <- cbind(g,l)
}

tab <- 
  replace_n %>%
  mutate(ID = 1:length(mg),
         sp = missing,
         sp.phy = NA) %>%
  arrange(genus)
tab

### Find species with potential replacement from the same genus on phy
c.mis <- 
  tab %>%
  filter(N > 0) %>%
  select(genus, sp)

## count missing sp per genus
n.miss.genus <- 
  c.mis %>%
  count(genus)

### Randomize new names based on species from the same genus--------------------
times <- 1000
sp.list <- list()
for (t in 1:times) {
  sp.phy <- list()
  for (k in 1:length(unique(tab$genus))) {
    i <- unique(tab$genus)[k]
    cc     <- filter(c.mis, genus == i)
    sp <- sample(replace[[i]], size = filter(n.miss.genus, genus == i)$n, replace = F)
    sp.phy[[k]] <- sp
  }
  sp.list[[t]] <- unlist(sp.phy)
}

### Final table with replaced names:
tab.list <- list()
for (i in 1:times) {
  tab$sp.phy <- sp.list[[i]]
  tab.list[[i]] <- tab
}
tab.list[[1]]

#check for duplicates with dataset species
dup <- numeric(times)
for(i in 1:times){
  dup[i] <- sum(as.character(d$sp) %in% tab.list[[i]]$sp.phy)
}

### Lists without duplicates
tab.list.ok <- tab.list[-which(dup > 0)]

### Randomize 300 from potential lists------------------------------------------
n_samp <- 300
samp <- sample(1:length(tab.list.ok), size = n_samp)
tab.list.ok <- tab.list.ok[samp]

### add species name to data
dmin.list <- list()
dminc <- d

for (i in 1:n_samp) {
  dminc$sp.phy <- as.character(dminc$sp)
  dminc[match(tab$sp, dminc$sp), ]$sp.phy <- tab.list.ok[[i]]$sp.phy
  rownames(dminc) <- as.character(dminc$sp.phy)
  dmin.list[[i]] <- dminc
}

### match data and phy (and replace names in the phylogenies)
d.list <- list()
for (i in 1:n_samp){
  d.list[[i]] <- match_dataphy(sp ~ 1, data = dmin.list[[i]], phy = p)
  d.list[[i]]$phy$tip.label  <- d.list[[i]]$data$sp
  rownames(d.list[[i]]$data) <- as.character(d.list[[i]]$data$sp)
}

### Test the number of trees are correct:
length(d.list) == n_samp

### Split data and tree:
trees <- list()
datas <- list()
for (i in 1:n_samp) {
  trees[[i]] <- d.list[[i]]$phy
  datas[[i]] <- d.list[[i]]$data
}
class(trees) <- "multiPhylo"

# Save the list of species replaced names on the phylogeny----------------------
w_repla <- which(datas[[1]]$sp != datas[[1]]$sp.phy)
tab_replaced_names <- data.frame(original_name = datas[[1]]$sp[w_repla])

for (i in seq_len(300)) {
  tab_replaced_names[, i + 1] <- datas[[i]][w_repla,]$sp.phy
}
names(tab_replaced_names) <- c("original_name", paste("replaced_on_tree_", 1:300,
                                                      "_by", sep = ""))
write_csv(x = tab_replaced_names, path = "outputs/temp/Table_replaced_names_for_missing_species.csv")

# Merge all trees and data
datas[[1]] <- datas[[1]] %>%  select(-sp.phy)
comp_data  <- match_dataphy(formula = tot ~ 1, phy = trees, data = datas[[1]])

# Randomize working tree (for representation and main analysis)-----------------
set.seed(12345)
n.tree <- sample(x = 1:300, size = 1)
wt <- trees[[n.tree]]

### Save 300 trees and data matched for analysis-----
write.tree(phy = trees, file = "data/processed/300_study_phylo_trees.tre")
write.tree(phy = wt, file = "data/processed/working_study_phylo_tree.tre")
write_csv(x = datas[[1]], path = "data/processed/data_flower_biomass_partition.csv")
saveRDS(object = comp_data, file = "data/processed/data_flower_biomass_with_tree.Rds")
### END-----
