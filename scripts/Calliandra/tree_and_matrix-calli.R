#=====================================================================================================#

library(ape)
library(tidyverse)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
calli_cr <- read.csv("datasets/Calliandra/calli_cr.csv", na.strings = c("", NA))

#Reading tree
calli_tree <- read.tree("trees/Calliandra/calli_tree-vasconcelos_etal.txt")

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree
calli_cr <- calli_cr[which(calli_cr$gen_sp %in% calli_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- calli_cr[ , which(colnames(calli_cr) %in% c("gen_sp", "id_grid"))]
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
calli_cr <- calli_cr[!calli_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

calli_pruned.tree <- drop.tip(calli_tree, 
                             calli_tree$tip.label[!calli_tree$tip.label %in% calli_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(calli_pruned.tree$edge.length)){
  edge.l[i] <- calli_pruned.tree$edge.length[i]/sum(calli_pruned.tree$edge.length)
}

calli_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(calli_pruned.tree, file = "trees/Calliandra/pruned_tree-calli.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
calli_matrix <- matrix(data = NA, nrow = length(unique(calli_cr$id_grid)), 
                         ncol = length(unique(calli_cr$gen_sp)))
calli_matrix <- as.data.frame(calli_matrix)
colnames(calli_matrix) <- unique(calli_cr$gen_sp)
rownames(calli_matrix) <- unique(calli_cr$id_grid)
for(i in 1:nrow(calli_matrix)){
  for(j in 1:ncol(calli_matrix)){
    if(colnames(calli_matrix)[j] %in% calli_cr$gen_sp[calli_cr$id_grid == rownames(calli_matrix)[i]]){
      calli_matrix[i, j] <- 1
    } else {
      calli_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(calli_matrix, "matrices/calli_matrix.csv", row.names = TRUE)
