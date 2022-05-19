#This code is intended to prune and manipulate the phylogenetic tree and create a 
#presence/absence matrix

#=====================================================================================================#

library(ape)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
vel_cr <- read.csv("datasets/Velloziaceae/vel_cr.csv", na.strings = c("", NA))

#Reading tree
vel_tree <- read.nexus("trees/Velloziaceae/velloziaceae_tree-vasconcelos_etal.nex")

#=====================================================================================================#

#====================#
# STANDARDIZING TREE #
#====================#

#Extracting labels
tip_labels <- data.frame("tip_labels" = vel_tree$tip.label, 
                         "std_labels" = NA, stringsAsFactors = FALSE)

#Standardizing tip labels
tip_labels[ , 2] <- trimws(tip_labels[ , 1])

#Replacing names in the tree
vel_tree$tip.label <- tip_labels$std_labels

#=================#
# COLLAPSING VARS #
#=================#

#Supressing infraspecific epithet
tip_labels$vars_sup <- NA
for(i in 1:nrow(tip_labels)){
  tip_labels[i, "vars_sup"] <- paste(strsplit(tip_labels[i, "std_labels"], "_")[[1]][c(1, 2)], 
                                     collapse = "_")
}
vel_tree$tip.label <- tip_labels$vars_sup

#Removing duplicated labels
duplicated_logical <- duplicated(vel_tree$tip.label)
vel_tree$tip.label <- as.character(1:length(vel_tree$tip.label))
for(i in 1:length(vel_tree$tip.label)){
  if(duplicated_logical[i] == TRUE)
    vel_tree <- drop.tip(vel_tree, vel_tree$tip.label[vel_tree$tip.label == as.character(i)])
}
vel_tree$tip.label <- tip_labels$vars_sup[-which(duplicated_logical)]

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree
vel_cr <- vel_cr[which(vel_cr$gen_sp %in% vel_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- vel_cr[ , which(colnames(vel_cr) %in% c("gen_sp", "id_grid"))]
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
vel_cr <- vel_cr[!vel_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

vel_pruned.tree <- drop.tip(vel_tree, 
                                   vel_tree$tip.label[!vel_tree$tip.label %in% vel_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(vel_pruned.tree$edge.length)){
  edge.l[i] <- vel_pruned.tree$edge.length[i]/sum(vel_pruned.tree$edge.length)
}

vel_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(vel_pruned.tree, file = "trees/Velloziaceae/pruned_tree-vel.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
vel_matrix <- matrix(data = NA, nrow = length(unique(vel_cr$id_grid)), 
                            ncol = length(unique(vel_cr$gen_sp)))
vel_matrix <- as.data.frame(vel_matrix)
colnames(vel_matrix) <- unique(vel_cr$gen_sp)
rownames(vel_matrix) <- unique(vel_cr$id_grid)
for(i in 1:nrow(vel_matrix)){
  for(j in 1:ncol(vel_matrix)){
    if(colnames(vel_matrix)[j] %in% vel_cr$gen_sp[vel_cr$id_grid == rownames(vel_matrix)[i]]){
      vel_matrix[i, j] <- 1
    } else {
      vel_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(vel_matrix, "matrices/vel_matrix.csv", row.names = TRUE)