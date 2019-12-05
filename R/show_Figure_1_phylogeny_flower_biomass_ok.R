# Figure 1: Phylogeny, maleness and Flower biomass

### Packages ----------------------------------------------------------------
library(sensiPhy)
library(phytools)
library(tidyverse)

### Start-----------------------------------------------------------------------
rm(list = ls())
source("R/ZZZ_functions.R")

### Load data-------------------------------------------------------------------
dall <- read.csv("data/processed/data_flower_biomass_partition.csv")
p    <- read.tree(file = "data/processed/working_study_phylo_tree.tre")

## Add species name as rows names:
row.names(dall) <- gsub(" ", "_", dall$sp)
dat1 <- sensiPhy::match_dataphy(tot ~ 1, data = dall, phy = p)

### Fig1A Biomass Histogram----
### Distribution of flower organs biomass
label      <- c("Male", "Female", "Petals", "Sepals")
colorblind <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")

### Plot parameters
bre    <- c(0.00001, 0.0001,0.001,0.01,0.1,1)
lab    <- c("0.00001", "0.0001","0.001","0.01","0.1","1")
lim    <- range(bre)
alp = .6

### Flower biomass:
fig1A <- 
  ggplot(dall, aes(x = tot)) + 
  geom_histogram(fill = gray(.2), alpha = alp, bins = 20,
                 position = "dodge", color = "black", size = .25) +
  scale_x_log10(breaks = bre, labels = lab) +
  xlab("Flower biomass (g)") +
  ylab("Number of Species") + 
  tema(base_size = 24);fig1A

### FLower biomass variation----
range(dall$tot)

# How many flower can be produced between the largest and the smallest flower?
range(dall$tot)[2]/range(dall$tot)[1]

### salvar:
ggsave("outputs/figures/raw/Fig_1A.pdf", 
       height = 6, width = 6, units = "in",
       plot = fig1A) 

### Fig1B Flower maleness----
dall$and.gyn <- log10(dall$and/dall$gyn)

colfunc <- colorRampPalette(c(colorblind[2], colorblind[1]))
cols1 <- colfunc(19)

fig1B <-
  ggplot(dall, aes(x = and.gyn, fill = and.gyn)) + 
  geom_histogram(alpha = .7, bins = 20, fill = cols1,
                 color = "black", size = .25) +
  geom_vline(xintercept = 0, lty = "dashed", color = "black") +
  scale_x_continuous(limits = c(-2, 2),
                     breaks = seq(-2, 2, 1),
                     labels = c("1:100", "1:10", "1:1", "10:1", "100:1")) +
  xlab("Male : Female ratio") +
  ylab("Number of species") + 
  tema(base_size = 24);fig1B

### salvar
ggsave("outputs/figures/raw/Fig_1B.pdf", 
       height = 6, width = 6, units = "in",
       plot = fig1B)

### Fig1C Phylogentic Tree----
X <- dat1$data
tree <- dat1$phy
X <- X %>%
  mutate(and.p = and/tot,
         gyn.p = gyn/tot,
         sep.p = sep/tot,
         pet.p = pet/tot) %>%
  select(sp, and.p, gyn.p, sep.p, pet.p)

rownames(X) <- X$sp
X <- select(X, -sp)
Y <- t(apply(X,1,cumsum))

cols  <- c("#D55E00","#56B4E9","#009E73", "#CC79A7")
scale <- 0.2 * max(nodeHeights(tree))
pdf(file = "outputs/figures/raw/Fig_1C.pdf")
plotTree.wBars(tree,Y[,ncol(Y)], type = "fan",
               scale = 50, col = cols[ncol(Y)],
               width = 3.5, border = "white", lwd = 1)

obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
for (i in (ncol(Y) - 1):1) {
  plotTree.wBars(tree, Y[,i], type = "fan", scale = 50, add = TRUE,
                 lims = obj$x.lim,
                 col = cols[i],
                 width = 3.5, 
                 border = "white",
                 lwd = .3,
                 color = gray(.1))
}
legend( x = "bottomright", legend = c("Male", "Female", "Sepals", "Petals"), pch = 22, 
        pt.cex = 2, pt.bg = cols, cex = 1.7, bty = "n", horiz = FALSE)
dev.off()
