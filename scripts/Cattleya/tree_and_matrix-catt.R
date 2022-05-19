#=====================================================================================================#

library(ape)
library(tidyverse)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
catt_cr <- read.csv("datasets/Cattleya/catt_cr.csv", na.strings = c("", NA))

#Reading tree
catt_tree <- read.tree("trees/Cattleya/catt_tree-vasconcelos_etal.txt")

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree
catt_cr <- catt_cr[which(catt_cr$gen_sp %in% catt_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- catt_cr[ , which(colnames(catt_cr) %in% c("gen_sp", "id_grid"))]
one_sample <- unique(all_samples, by = c("gen_sp", "id_grid"))

#Redundancy values
redundancy <- tibble("id_grid" = unique(all_samples$id_grid), "richness" = NA, "n" = NA,
                     "redundancy" = NA)

SR <- plyr::count(one_sample$id_grid)
N <- plyr::count(all_samples$id_grid)

for(i in 1:nrow(redundancy)){
  redundancy$richness[i] <- SR$freq[SR$x == redundancy$id_grid[i]]
  redundancy$n[i] <- N$freq[N$x == redundancy$id_grid[i]]
}

redundancy$redundancy <- 1-(redundancy$richness/redundancy$n)

#Removing grids with redundancy = 0 (meaning a species per sample)
inv_grids <- redundancy$id_grid[redundancy$redundancy == 0]
catt_cr <- catt_cr[!catt_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

catt_pruned.tree <- drop.tip(catt_tree, 
                              catt_tree$tip.label[!catt_tree$tip.label %in% catt_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(catt_pruned.tree$edge.length)){
  edge.l[i] <- catt_pruned.tree$edge.length[i]/sum(catt_pruned.tree$edge.length)
}

catt_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(catt_pruned.tree, file = "trees/Cattleya/pruned_tree-catt.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
catt_matrix <- matrix(data = NA, nrow = length(unique(catt_cr$id_grid)), 
                       ncol = length(unique(catt_cr$gen_sp)))
catt_matrix <- as.data.frame(catt_matrix)
colnames(catt_matrix) <- unique(catt_cr$gen_sp)
rownames(catt_matrix) <- unique(catt_cr$id_grid)
for(i in 1:nrow(catt_matrix)){
  for(j in 1:ncol(catt_matrix)){
    if(colnames(catt_matrix)[j] %in% catt_cr$gen_sp[catt_cr$id_grid == rownames(catt_matrix)[i]]){
      catt_matrix[i, j] <- 1
    } else {
      catt_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(catt_matrix, "matrices/catt_matrix.csv", row.names = TRUE)
