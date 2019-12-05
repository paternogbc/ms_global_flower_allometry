# Supplementary material
# Figure: Flower and flower organs (distribution of biomass)

# Packages ----------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(cowplot)
source(file = "R/ZZZ_functions.R")

# Load data ---------------------------------------------------------------
d <- read.csv("data/processed/data_flower_biomass_partition.csv")

label      <- c("Male", "Female", "Petals", "Sepals")
colorblind <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")

### Plot parameters
bre    <- c(0.00001, 0.0001,0.001,0.01,0.1,1)
lab    <- c("0.00001", "0.0001","0.001","0.01","0.1","1")
lim    <- range(bre)
alp = .6
axi.tex = 14
bsize = 10

### Flower
ball <- 
  ggplot(d, aes(x = tot)) + 
  geom_histogram(fill = "black", alpha = alp, bins = 30,
                 position = "dodge", color = "black", size = .25) +
  scale_x_log10(breaks = bre, labels = lab) +
  xlab("Flower biomass (g)") +
  ylab("Frequency") + 
  theme_bw(base_size = bsize ) + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = axi.tex));ball

ball <-
  ggdraw() +
  draw_plot(ball) +
  draw_image(image = "images/flower_ilustrations/flower.png", 
             x = .35, 
             y = .32, 
             scale = .25)

  ### Male biomass
band <- 
  ggplot(d, aes(x = and)) + 
  geom_histogram(fill = colorblind[1], alpha = alp, 
                 position = "dodge", bins = 30, color = colorblind[1],
                 size = .25) +
  scale_x_log10(breaks = bre, labels = lab) +
  xlab("Male biomass (g)") +
  ylab("Frequency") + 
  theme_bw(base_size = bsize ) + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = axi.tex));band

band <-
  ggdraw() +
  draw_plot(band) +
  draw_image(image = "images/flower_ilustrations/male.png", 
             x = .35, 
             y = .32, 
             scale = .25);band

### Female biomass
bgyn <- 
  ggplot(d, aes(x = gyn)) + 
  geom_histogram(fill = colorblind[2], alpha = alp,
                 position = "dodge", bins = 30,
                 color = colorblind[2],
                 size = .25) +
  scale_x_log10(breaks = bre, labels = lab) +
  xlab("Female biomass (g)") +
  ylab("Frequency") + 
  theme_bw(base_size = bsize ) + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = axi.tex));bgyn

bgyn <-
  ggdraw() +
  draw_plot(bgyn) +
  draw_image(image = "images/flower_ilustrations/female.png", 
             x = .35, 
             y = .32, 
             scale = .25);bgyn
### Petals biomass
bpet <- 
  ggplot(d, aes(x = pet)) + 
  geom_histogram(fill = colorblind[3], alpha = alp,
                 position = "dodge", bins = 30,
                 color = colorblind[3],
                 size = .25) +
  scale_x_log10(breaks = bre, labels = lab) +
  xlab("Petals biomass (g)") +
  ylab("Frequency") + 
  theme_bw(base_size = bsize ) + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = axi.tex));bpet
bpet <-
  ggdraw() +
  draw_plot(bpet) +
  draw_image(image = "images/flower_ilustrations/petals.png", 
             x = .35, 
             y = .32, 
             scale = .25);bpet

### Sepals biomass
bsep <- 
  ggplot(d, aes(x = sep)) + 
  geom_histogram(fill = colorblind[4], alpha = alp,
                 position = "dodge", bins = 30,
                 color = colorblind[4],
                 size = .25) +
  scale_x_log10(breaks = bre, labels = lab) +
  xlab("Sepals biomass (g)") +
  ylab("Frequency") + 
  theme_bw(base_size = bsize ) + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = axi.tex));bsep

bsep <-
  ggdraw() +
  draw_plot(bsep) +
  draw_image(image = "images/flower_ilustrations/sepals.png", 
             x = .35, 
             y = .32, 
             scale = .25);bsep

### Arrange figure S3
g1 <- 
  plot_grid(ball, band, bgyn, bpet, bsep,  
            labels = c("A", "B", "C", "D", "E"), label_fontface = "plain",
            ncol = 1, nrow = 5);g1

### Export Figure
ggsave("outputs/figures/supp/SFig_flower_biomass_distribution.png", 
       height = 12, width = 4, units = "in",
       plot = g1)
