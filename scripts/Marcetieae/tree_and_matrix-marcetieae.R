#=====================================================================================================#

library(ape)
library(tidyverse)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
marcetieae_cr <- read.csv("datasets/Marcetieae/marcetieae_cr.csv", na.strings = c("", NA))

#Reading tree
marcetieae_tree <- read.tree("trees/Marcetieae/marcetieae_tree-vasconcelos_etal.txt")

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree (2,438)
marcetieae_cr <- marcetieae_cr[which(marcetieae_cr$gen_sp %in% marcetieae_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- marcetieae_cr[ , which(colnames(marcetieae_cr) %in% c("gen_sp", "id_grid"))]
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

#Removing grids with redundancy = 0 (meaning a species per sample; 2,421)
inv_grids <- redundancy$id_grid[redundancy$redundancy == 0]
marcetieae_cr <- marcetieae_cr[!marcetieae_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

marcetieae_pruned.tree <- drop.tip(marcetieae_tree, 
                                   marcetieae_tree$tip.label[!marcetieae_tree$tip.label %in% marcetieae_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(marcetieae_pruned.tree$edge.length)){
  edge.l[i] <- marcetieae_pruned.tree$edge.length[i]/sum(marcetieae_pruned.tree$edge.length)
}

marcetieae_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(marcetieae_pruned.tree, file = "trees/Marcetieae/pruned_tree-marcetieae.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
marcetieae_matrix <- matrix(data = NA, nrow = length(unique(marcetieae_cr$id_grid)), 
                            ncol = length(unique(marcetieae_cr$gen_sp)))
marcetieae_matrix <- as.data.frame(marcetieae_matrix)
colnames(marcetieae_matrix) <- unique(marcetieae_cr$gen_sp)
rownames(marcetieae_matrix) <- unique(marcetieae_cr$id_grid)
for(i in 1:nrow(marcetieae_matrix)){
  for(j in 1:ncol(marcetieae_matrix)){
    if(colnames(marcetieae_matrix)[j] %in% marcetieae_cr$gen_sp[marcetieae_cr$id_grid == rownames(marcetieae_matrix)[i]]){
      marcetieae_matrix[i, j] <- 1
    } else {
      marcetieae_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(marcetieae_matrix, "matrices/marcetieae_matrix.csv", row.names = TRUE)
