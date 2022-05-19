#=====================================================================================================#

library(ape)
library(tidyverse)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
minaria_cr <- read.csv("datasets/Minaria/minaria_cr.csv", na.strings = c("", NA))

#Reading tree
minaria_tree <- read.tree("trees/Minaria/minaria_tree-vasconcelos_etal.txt")

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree
minaria_cr <- minaria_cr[which(minaria_cr$gen_sp %in% minaria_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- minaria_cr[ , which(colnames(minaria_cr) %in% c("gen_sp", "id_grid"))]
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
minaria_cr <- minaria_cr[!minaria_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

minaria_pruned.tree <- drop.tip(minaria_tree, 
                              minaria_tree$tip.label[!minaria_tree$tip.label %in% minaria_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(minaria_pruned.tree$edge.length)){
  edge.l[i] <- minaria_pruned.tree$edge.length[i]/sum(minaria_pruned.tree$edge.length)
}

minaria_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(minaria_pruned.tree, file = "trees/Minaria/pruned_tree-minaria.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
minaria_matrix <- matrix(data = NA, nrow = length(unique(minaria_cr$id_grid)), 
                        ncol = length(unique(minaria_cr$gen_sp)))
minaria_matrix <- as.data.frame(minaria_matrix)
colnames(minaria_matrix) <- unique(minaria_cr$gen_sp)
rownames(minaria_matrix) <- unique(minaria_cr$id_grid)
for(i in 1:nrow(minaria_matrix)){
  for(j in 1:ncol(minaria_matrix)){
    if(colnames(minaria_matrix)[j] %in% minaria_cr$gen_sp[minaria_cr$id_grid == rownames(minaria_matrix)[i]]){
      minaria_matrix[i, j] <- 1
    } else {
      minaria_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(minaria_matrix, "matrices/minaria_matrix.csv", row.names = TRUE)
