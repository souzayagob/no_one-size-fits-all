#This code is intended to prune and manipulate the phylogenetic tree and create a 
#presence/absence matrix

#=====================================================================================================#

library(ape)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
diplusodon_cr <- read.csv("datasets/Diplusodon/diplusodon_cr.csv", na.strings = c("", NA))

#Reading tree
diplusodon_tree <- read.nexus("trees/Diplusodon/diplusodon_tree-vasconcelos_etal.nex")

#=====================================================================================================#

#====================#
# STANDARDIZING TREE #
#====================#

#Extracting tip labels
tip_labels <- data.frame("tip_labels" = diplusodon_tree$tip.label, 
                         "std_labels" = NA, stringsAsFactors = FALSE)

#Standardizing tip labels
for(i in 1:nrow(tip_labels)){
  tip_labels[i, 2] <- gsub("var_", "", tip_labels[i, 1])
  tip_labels[i, 2] <- trimws(tip_labels[i, 2])
} 

#Replacing names in the tree
diplusodon_tree$tip.label <- tip_labels$std_labels

#=================#
# COLLAPSING VARS #
#=================#

#Supressing infraspecific epithet
tip_labels$vars_sup <- NA
for(i in 1:nrow(tip_labels)){
  tip_labels[i, "vars_sup"] <- paste(strsplit(tip_labels[i, "std_labels"], "_")[[1]][c(1, 2)], 
                                     collapse = "_")
}
diplusodon_tree$tip.label <- tip_labels$vars_sup

#Removing duplicated labels
duplicated_logical <- duplicated(diplusodon_tree$tip.label)
diplusodon_tree$tip.label <- as.character(1:length(diplusodon_tree$tip.label))
for(i in 1:length(diplusodon_tree$tip.label)){
  if(duplicated_logical[i] == TRUE)
    diplusodon_tree <- drop.tip(diplusodon_tree, diplusodon_tree$tip.label[diplusodon_tree$tip.label == as.character(i)])
}
diplusodon_tree$tip.label <- tip_labels$vars_sup[-which(duplicated_logical)]

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree
diplusodon_cr <- diplusodon_cr[which(diplusodon_cr$gen_sp %in% diplusodon_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- diplusodon_cr[ , which(colnames(diplusodon_cr) %in% c("gen_sp", "id_grid"))]
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
diplusodon_cr <- diplusodon_cr[!diplusodon_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

diplusodon_pruned.tree <- drop.tip(diplusodon_tree, 
                                   diplusodon_tree$tip.label[!diplusodon_tree$tip.label %in% diplusodon_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(diplusodon_pruned.tree$edge.length)){
  edge.l[i] <- diplusodon_pruned.tree$edge.length[i]/sum(diplusodon_pruned.tree$edge.length)
}

diplusodon_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(diplusodon_pruned.tree, file = "trees/Diplusodon/pruned_tree-diplusodon.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
diplusodon_matrix <- matrix(data = NA, nrow = length(unique(diplusodon_cr$id_grid)), 
                        ncol = length(unique(diplusodon_cr$gen_sp)))
diplusodon_matrix <- as.data.frame(diplusodon_matrix)
colnames(diplusodon_matrix) <- unique(diplusodon_cr$gen_sp)
rownames(diplusodon_matrix) <- unique(diplusodon_cr$id_grid)
for(i in 1:nrow(diplusodon_matrix)){
  for(j in 1:ncol(diplusodon_matrix)){
    if(colnames(diplusodon_matrix)[j] %in% diplusodon_cr$gen_sp[diplusodon_cr$id_grid == rownames(diplusodon_matrix)[i]]){
      diplusodon_matrix[i, j] <- 1
    } else {
      diplusodon_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(diplusodon_matrix, "matrices/diplusodon_matrix.csv", row.names = TRUE)
