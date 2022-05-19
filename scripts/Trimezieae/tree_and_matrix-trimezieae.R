#=====================================================================================================#

library(ape)
library(tidyverse)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
trimezieae_cr <- read.csv("datasets/Trimezieae/trimezieae_cr.csv", na.strings = c("", NA))

#Reading tree
trimezieae_tree <- read.tree("trees/Trimezieae/trimezieae_tree-vasconcelos_etal.txt")

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree (496)
trimezieae_cr <- trimezieae_cr[which(trimezieae_cr$gen_sp %in% trimezieae_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- trimezieae_cr[ , which(colnames(trimezieae_cr) %in% c("gen_sp", "id_grid"))]
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

#Removing grids with redundancy = 0 (meaning a species per sample; 491)
inv_grids <- redundancy$id_grid[redundancy$redundancy == 0]
trimezieae_cr <- trimezieae_cr[!trimezieae_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

trimezieae_pruned.tree <- drop.tip(trimezieae_tree, 
                                  trimezieae_tree$tip.label[!trimezieae_tree$tip.label %in% trimezieae_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(trimezieae_pruned.tree$edge.length)){
  edge.l[i] <- trimezieae_pruned.tree$edge.length[i]/sum(trimezieae_pruned.tree$edge.length)
}

trimezieae_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(trimezieae_pruned.tree, file = "trees/Trimezieae/pruned_tree-trimezieae.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
trimezieae_matrix <- matrix(data = NA, nrow = length(unique(trimezieae_cr$id_grid)), 
                           ncol = length(unique(trimezieae_cr$gen_sp)))
trimezieae_matrix <- as.data.frame(trimezieae_matrix)
colnames(trimezieae_matrix) <- unique(trimezieae_cr$gen_sp)
rownames(trimezieae_matrix) <- unique(trimezieae_cr$id_grid)
for(i in 1:nrow(trimezieae_matrix)){
  for(j in 1:ncol(trimezieae_matrix)){
    if(colnames(trimezieae_matrix)[j] %in% trimezieae_cr$gen_sp[trimezieae_cr$id_grid == rownames(trimezieae_matrix)[i]]){
      trimezieae_matrix[i, j] <- 1
    } else {
      trimezieae_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(trimezieae_matrix, "matrices/trimezieae_matrix.csv", row.names = TRUE)
