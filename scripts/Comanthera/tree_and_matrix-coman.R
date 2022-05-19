#=====================================================================================================#

library(ape)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
coman_cr <- read.csv("datasets/Comanthera/coman_cr.csv", na.strings = c("", NA),
                     stringsAsFactors = FALSE)

#Reading tree
coman_tree <- read.nexus("trees/Comanthera/comanthera_tree.tre")[[1]]

#=====================================================================================================#

#====================#
# STANDARDIZING TREE #
#====================#

#Extracting labels
tip_labels <- data.frame("tip_labels" = coman_tree$tip.label, 
                         "std_labels" = NA, "corrected_labels" = NA)

#Standardizing labels
for(i in 1:nrow(tip_labels)){
  tip_labels[i, 2] <- gsub("C", "Comanthera_", tip_labels[i, 1])
  tip_labels[i, 2] <- trimws(tip_labels[i, 2])
} 

#Replacing names in the tree
coman_tree$tip.label <- tip_labels$std_labels

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree
coman_cr <- coman_cr[which(coman_cr$gen_sp %in% coman_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- coman_cr[ , which(colnames(coman_cr) %in% c("gen_sp", "id_grid"))]
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
coman_cr <- coman_cr[!coman_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

coman_pruned.tree <- drop.tip(coman_tree, 
                                   coman_tree$tip.label[!coman_tree$tip.label %in% coman_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(coman_pruned.tree$edge.length)){
  edge.l[i] <- coman_pruned.tree$edge.length[i]/sum(coman_pruned.tree$edge.length)
}

coman_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

#write.nexus(coman_pruned.tree, file = "trees/Comanthera/pruned_tree-coman.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
coman_matrix <- matrix(data = NA, nrow = length(unique(coman_cr$id_grid)), 
                        ncol = length(unique(coman_cr$gen_sp)))
coman_matrix <- as.data.frame(coman_matrix)
colnames(coman_matrix) <- unique(coman_cr$gen_sp)
rownames(coman_matrix) <- unique(coman_cr$id_grid)
for(i in 1:nrow(coman_matrix)){
  for(j in 1:ncol(coman_matrix)){
    if(colnames(coman_matrix)[j] %in% coman_cr$gen_sp[coman_cr$id_grid == rownames(coman_matrix)[i]]){
      coman_matrix[i, j] <- 1
    } else {
      coman_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(coman_matrix, "matrices/coman_matrix.csv", row.names = TRUE)
