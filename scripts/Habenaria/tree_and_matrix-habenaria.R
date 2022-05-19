#=====================================================================================================#

library(ape)
library(tidyverse)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
habenaria_cr <- read.csv("datasets/Habenaria/habenaria_cr.csv", na.strings = c("", NA))

#Reading tree
habenaria_tree <- read.tree("trees/Habenaria/habenaria_tree-vasconcelos_etal.txt")

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree (780)
habenaria_cr <- habenaria_cr[which(habenaria_cr$gen_sp %in% habenaria_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- habenaria_cr[ , which(colnames(habenaria_cr) %in% c("gen_sp", "id_grid"))]
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

#Removing grids with redundancy = 0 (meaning a species per sample; 773)
inv_grids <- redundancy$id_grid[redundancy$redundancy == 0]
habenaria_cr <- habenaria_cr[!habenaria_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

habenaria_pruned.tree <- drop.tip(habenaria_tree, 
                             habenaria_tree$tip.label[!habenaria_tree$tip.label %in% habenaria_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(habenaria_pruned.tree$edge.length)){
  edge.l[i] <- habenaria_pruned.tree$edge.length[i]/sum(habenaria_pruned.tree$edge.length)
}

habenaria_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(habenaria_pruned.tree, file = "trees/Habenaria/pruned_tree-habenaria.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
habenaria_matrix <- matrix(data = NA, nrow = length(unique(habenaria_cr$id_grid)), 
                      ncol = length(unique(habenaria_cr$gen_sp)))
habenaria_matrix <- as.data.frame(habenaria_matrix)
colnames(habenaria_matrix) <- unique(habenaria_cr$gen_sp)
rownames(habenaria_matrix) <- unique(habenaria_cr$id_grid)
for(i in 1:nrow(habenaria_matrix)){
  for(j in 1:ncol(habenaria_matrix)){
    if(colnames(habenaria_matrix)[j] %in% habenaria_cr$gen_sp[habenaria_cr$id_grid == rownames(habenaria_matrix)[i]]){
      habenaria_matrix[i, j] <- 1
    } else {
      habenaria_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(habenaria_matrix, "matrices/habenaria_matrix.csv", row.names = TRUE)
