#This code is intended to prune and manipulate the phylogenetic tree and create a 
#presence/absence matrix

#=====================================================================================================#

library(ape)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
paep_cr <- read.csv("datasets/Paepalanthus/paep_cr.csv", na.strings = c("", NA))

#Reading tree
paep_tree <- read.nexus("trees/Paepalanthus/Paepalanthus_Bayes.nex")

#=====================================================================================================#

#====================#
# STANDARDIZING TREE #
#====================#

#Extracting tip labels
tip_labels <- data.frame("tip_labels" = paep_tree$tip.label, 
                         "std_labels" = NA, stringsAsFactors = FALSE)

#Standardizing tip labels (maintaining varieties)
vars <- tip_labels[grep("var_", tip_labels$tip_labels),  1] 
for(i in 1:nrow(tip_labels)){
  if(!tip_labels$tip_labels[i] %in% vars){
    tip_labels[i, 2] <- paste(strsplit(tip_labels[i, 1], 
                                       split = "_")[[1]][c(1,2)], collapse = "_")
  } else{
    tip_labels[i, 2] <- paste(strsplit(tip_labels[i, 1], 
                                       split = "_")[[1]][c(1,2,3,4)], collapse = "_")
  }
}
tip_labels$std_labels <- gsub("_var_", "_", tip_labels$std_labels)

#Replacing tip labels in the tree
paep_tree$tip.label <- tip_labels$std_labels

#=================#
# COLLAPSING VARS #
#=================#

#Supressing infraspecific epithet
tip_labels$vars_sup <- NA
for(i in 1:nrow(tip_labels)){
  tip_labels[i, "vars_sup"] <- paste(strsplit(tip_labels[i, "std_labels"], "_")[[1]][c(1, 2)], 
                                     collapse = "_")
}
paep_tree$tip.label <- tip_labels$vars_sup

#Removing duplicated labels
duplicated_logical <- duplicated(paep_tree$tip.label)
paep_tree$tip.label <- as.character(1:length(paep_tree$tip.label))
for(i in 1:length(paep_tree$tip.label)){
  if(duplicated_logical[i] == TRUE)
    paep_tree <- drop.tip(paep_tree, paep_tree$tip.label[paep_tree$tip.label == as.character(i)])
}
paep_tree$tip.label <- tip_labels$vars_sup[-which(duplicated_logical)]

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree
paep_cr <- paep_cr[which(paep_cr$gen_sp %in% paep_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- paep_cr[ , which(colnames(paep_cr) %in% c("gen_sp", "id_grid"))]
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
paep_cr <- paep_cr[!paep_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

paep_pruned.tree <- drop.tip(paep_tree, 
                                   paep_tree$tip.label[!paep_tree$tip.label %in% paep_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(paep_pruned.tree$edge.length)){
  edge.l[i] <- paep_pruned.tree$edge.length[i]/sum(paep_pruned.tree$edge.length)
}

paep_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(paep_pruned.tree, file = "trees/Paepalanthus/pruned_tree-paep.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
paep_matrix <- matrix(data = NA, nrow = length(unique(paep_cr$id_grid)), 
                            ncol = length(unique(paep_cr$gen_sp)))
paep_matrix <- as.data.frame(paep_matrix)
colnames(paep_matrix) <- unique(paep_cr$gen_sp)
rownames(paep_matrix) <- unique(paep_cr$id_grid)
for(i in 1:nrow(paep_matrix)){
  for(j in 1:ncol(paep_matrix)){
    if(colnames(paep_matrix)[j] %in% paep_cr$gen_sp[paep_cr$id_grid == rownames(paep_matrix)[i]]){
      paep_matrix[i, j] <- 1
    } else {
      paep_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(paep_matrix, "matrices/paep_matrix.csv", row.names = TRUE)
