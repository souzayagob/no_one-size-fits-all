#=====================================================================================================#

library(raster)
library(rgdal)
library(tidyverse)
library(ape)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

#=====================================================================================================#

#=======#
# INPUT #
#=======#

# Reading clean dataset
coords <- read.csv(file = "datasets/Minaria/minaria_cr.csv", na.strings = c("", NA))

# Reading PD results
pd <- read.csv("results/Minaria/pd_stats.csv")
pd <- pd %>% rename("id_grid" = id)

# Reading tree
tree <- read.nexus("trees/Minaria/pruned_tree-minaria.nex")

# Removing from the dataset all species that are not represented in the tree
coords <- coords[coords$gen_sp %in% tree$tip.label, ]

# Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- coords[ , which(colnames(coords) %in% c("gen_sp", "id_grid"))]
one_sample <- unique(all_samples, by = c("gen_sp", "id_grid"))

# Redundancy values
redundancy <- tibble("id_grid" = unique(all_samples$id_grid), "richness" = NA, "n" = NA,
                     "redundancy" = NA)

SR <- plyr::count(one_sample$id_grid)
N <- plyr::count(all_samples$id_grid)

for(i in 1:nrow(redundancy)){
  redundancy$richness[i] <- SR$freq[SR$x == redundancy$id_grid[i]]
  redundancy$n[i] <- N$freq[N$x == redundancy$id_grid[i]]
}

redundancy$redundancy <- 1-(redundancy$richness/redundancy$n)

# Merging redundancy and PD
pd_red <- merge(pd, redundancy, by = "id_grid")

corrPdRed <- as.character(formatC(cor(pd_red$pd.obs, 
                                      pd_red$redundancy, use = "na.or.complete"))) 
corrPdRed_plot <- ggplot(data = pd_red, mapping = aes(jitter(redundancy), pd.obs))+
  geom_jitter()+
  labs(title = "Minaria", subtitle = paste("r =", corrPdRed))+
  xlab("Redundancy")+
  ylab("Phylogenetic diversity")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))
