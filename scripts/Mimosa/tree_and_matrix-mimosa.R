#This code is intended to prune and manipulate the phylogenetic tree and create a 
#presence/absence matrix

#OBS.:
#(1) Here, I supress infraspecif epithets. However, it is more appropriate, I think, to include them.
#In order to this, I should review the code I used to clean taxonomic names

#=====================================================================================================#

library(ape)
library(tidyverse)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
mimosa_cr <- read.csv("datasets/Mimosa/mimosa_cr.csv", na.strings = c("", NA))

#Reading tree
mimosa_tree <- read.tree("trees/Mimosa/mimosa_tree-vasconcelos_etal.txt")

#=====================================================================================================#

#====================#
# STANDARDIZING TREE #
#====================#

#Extracting tip labels
tip_labels <- data.frame("tip_labels" =  mimosa_tree$tip.label, 
                         "std_labels" = NA)

#Standardizing tip labels
tip_labels$std_labels <- trimws(tip_labels$tip_labels)

#Replacing names in the tree
mimosa_tree$tip.label <- tip_labels$std_labels

#=================#
# COLLAPSING VARS #
#=================#

#Supressing infraspecific epithet
tip_labels$vars_sup <- NA
for(i in 1:nrow(tip_labels)){
  tip_labels[i, "vars_sup"] <- paste(strsplit(tip_labels[i, "std_labels"], "_")[[1]][c(1, 2)], 
                                     collapse = "_")
}
mimosa_tree$tip.label <- tip_labels$vars_sup

#Removing duplicated labels
duplicated_logical <- duplicated(mimosa_tree$tip.label)
mimosa_tree$tip.label <- as.character(1:length(mimosa_tree$tip.label))
for(i in 1:length(mimosa_tree$tip.label)){
  if(duplicated_logical[i] == TRUE)
  mimosa_tree <- drop.tip(mimosa_tree, mimosa_tree$tip.label[mimosa_tree$tip.label == as.character(i)])
}
mimosa_tree$tip.label <- tip_labels$vars_sup[-which(duplicated_logical)]

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree
mimosa_cr <- mimosa_cr[which(mimosa_cr$gen_sp %in% mimosa_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- mimosa_cr[ , which(colnames(mimosa_cr) %in% c("gen_sp", "id_grid"))]
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
mimosa_cr <- mimosa_cr[!mimosa_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

mimosa_pruned.tree <- drop.tip(mimosa_tree, mimosa_tree$tip.label[!mimosa_tree$tip.label %in% mimosa_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

#Pruned tree
edge.l <- c()
for(i in 1:length(mimosa_pruned.tree$edge.length)){
  edge.l[i] <- mimosa_pruned.tree$edge.length[i]/sum(mimosa_pruned.tree$edge.length)
}

mimosa_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(mimosa_pruned.tree, file = "trees/Mimosa/pruned_tree-mimosa.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
mimosa_matrix <- matrix(data = NA, nrow = length(unique(mimosa_cr$id_grid)), 
                        ncol = length(unique(mimosa_cr$gen_sp)))
mimosa_matrix <- as.data.frame(mimosa_matrix)
colnames(mimosa_matrix) <- unique(mimosa_cr$gen_sp)
rownames(mimosa_matrix) <- unique(mimosa_cr$id_grid)
for(i in 1:nrow(mimosa_matrix)){
  for(j in 1:ncol(mimosa_matrix)){
    if(colnames(mimosa_matrix)[j] %in% mimosa_cr$gen_sp[mimosa_cr$id_grid == rownames(mimosa_matrix)[i]]){
      mimosa_matrix[i, j] <- 1
    } else {
      mimosa_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(mimosa_matrix, "matrices/mimosa_matrix.csv", row.names = TRUE)
