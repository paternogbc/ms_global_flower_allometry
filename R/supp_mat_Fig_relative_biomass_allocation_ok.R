# Supplementary material
# Figure: Relative biomass allocation

# Packages ----------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(reshape2)
source(file = "R/ZZZ_functions.R")

# Load data ---------------------------------------------------------------
d <- read.csv("data/processed/data_flower_biomass_partition.csv")

raw.rel <- mutate(d, 
                  and.p = (and/tot)*100,
                  gyn.p = (gyn/tot)*100,
                  pet.p = (pet/tot)*100,
                  sep.p = (sep/tot)*100) %>% 
  select(and.p, gyn.p, pet.p, sep.p)

raw.rel.long <- melt(data = raw.rel, 
                     variable.name = "component",value.name = "biomass" )

### Allocated biomass across flower functions:
### Define color blind colors:
cols <- c("#D55E00","#56B4E9","#CC79A7", "#009E73")

grelative <- 
  ggplot(raw.rel.long, aes(y = biomass, x = component, fill = component)) + 
  geom_jitter(size = .5, width = .05, show.legend = F, color = gray(.5)) +
  geom_boxplot(show.legend = F, alpha = .75, 
               outlier.color = NA, color = gray(.2)) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  scale_x_discrete(labels = c("Male", "Female", "Petals", "Sepals")) +
  labs(y = "Percentage of biomass (%)", x = "Flower component") + 
  tema(base_size = 22);grelative

grelative <-
ggdraw() +
  draw_plot(grelative) +
  draw_image(
    image = "images/flower_ilustrations/male.png",
    x = -.22,
    y = .4,
    scale = .13
  ) +
  draw_image(
    image = "images/flower_ilustrations/female.png",
    x = -0.03,
    y = .4,
    scale = .13
  ) +
  draw_image(
    image = "images/flower_ilustrations/petals.png",
    x = .17,
    y = .4,
    scale = .13
  ) +
  draw_image(
    image = "images/flower_ilustrations/sepals.png",
    x = .36,
    y = .4,
    scale = .13
  )
grelative

### Figure ----
ggsave("outputs/figures/supp/SFig_relative_biomass_allocation.png", 
       height = 8, width = 6, units = "in",
       plot = grelative)

### Table with relative allocation summary
tab.rel.sum <-
  raw.rel.long %>%
  group_by(component) %>%
  summarise("Mean (%)" = mean(biomass),
            sd = sd(biomass),
            min = min(biomass),
            max = max(biomass)) %>% 
  rename("Flower component" = component)
tab.rel.sum
write_csv(tab.rel.sum, "outputs/tables/supp/STable_relative_biomass_allocation.csv")
### END----
