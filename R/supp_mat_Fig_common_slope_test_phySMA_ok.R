# Supplementary material
# Figure: common slope test

# Packages ----------------------------------------------------------------
library(tidyverse)
library(cowplot)
source(file = "R/ZZZ_functions.R")

# Load data ---------------------------------------------------------------
# Permutation test for phySMA
d <- readRDS("outputs/temp/common_slope_test_cache_phySMA.Rds")

# Male vs Female-----
# phySMA
g1 <- 
  ggplot(data = d[[1]]$null_distribution, aes(x = diff_null)) +
  geom_histogram(fill = "gray", color = "black") +
  geom_vline(xintercept = d[[1]]$slope_test$diff, color = "red") +
  tema() +
  labs(y = "Frequency", x = "Slope difference (null distribution)",
       title = "Common slope permutation test (phySMA)",
       subtitle = "Male vs Female",
       caption = paste("N of sim. = ", d[[1]]$slope_test$n_sim, 
                       "|","P.value = ", d[[1]]$slope_test$p.value))
g1

# petals vs Sepals-----
# phySMA
g2 <- 
  ggplot(data = d[[2]]$null_distribution, aes(x = diff_null)) +
  geom_histogram(fill = "gray", color = "black") +
  geom_vline(xintercept = d[[2]]$slope_test$diff, color = "red") +
  tema() +
  labs(y = "Frequency", x = "Slope difference (null distribution)",
       title = "Common slope permutation test (phySMA)",
       subtitle = "Petals vs Sepals",
       caption = paste("N of sim. = ", d[[2]]$slope_test$n_sim, 
                       "|","P.value = ", d[[2]]$slope_test$p.value))
g2

gall <- plot_grid(g1, g2, labels = c("A", "B"), label_size = 18, label_fontface = "plain")
gall

# Save plot----------
ggsave(plot = gall, filename = "outputs/figures/supp/SFig_common_slope_test_phySMA.png",
      width = 11, height = 5)
# END----