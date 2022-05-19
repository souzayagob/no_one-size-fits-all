#This code is intended to prune and manipulate the phylogenetic tree and create a 
#presence/absence matrix

#=====================================================================================================#

library(ape)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading CR dataset
lych_cr <- read.csv("datasets/Lychnophorinae/lych_cr.csv", na.strings = c("", NA))

#Reading tree
lych_tree <- read.nexus("trees/Lychnophorinae/lychnophorinae_tree-vasconcelos_etal.nex")

#=====================================================================================================#

#====================#
# STANDARDIZING TREE #
#====================#

#Extracting tip labels
tip_labels <- data.frame("tip_labels" = lych_tree$tip.label, 
                         "std_labels" = NA, stringsAsFactors = FALSE)

#Standardizing tip labels
for(i in 1:nrow(tip_labels)){
  tip_labels[i, 2] <- gsub("var_", "", tip_labels[i, 1])
  tip_labels[i, 2] <- trimws(tip_labels[i, 2])
} 

#Replacing names in the tree
lych_tree$tip.label <- tip_labels$std_labels

#=====================================================================================================#

#===============#
# SAMPLING BIAS #
#===============#

#Removing species not sampled in the tree
lych_cr <- lych_cr[which(lych_cr$gen_sp %in% lych_tree$tip.label), ]

#============#
# Redundancy #
#============#

#Data frames with all samples (all_samples) and with one sample per species per grid (one_sample)
all_samples <- lych_cr[ , which(colnames(lych_cr) %in% c("gen_sp", "id_grid"))]
one_sample <- unique(all_samples, by = c("gen_sp", "id_grid"))

#Redundancy values
redundancy <- tibble("id_grid" = unique(all_samples$id_grid), "richness" = NA, "n" = NA,
                     "redundancy" = NA)

SR <- plyr::count(one_sample$id_grid)
N <- plyr::count(all_samples$id_grid)

for(i in 1:nrow(redundancy)){
  redundancy$richness[i] <- SR$freq[SR$x == redunmellobarretoidancy$id_grid[i]]
  redundancy$n[i] <- N$freq[N$x == redundancy$id_grid[i]]
}

redundancy$redundancy <- 1-(redundancy$richness/redundancy$n)

#Removing grids with redundancy = 0 (meaning a species per sample)
inv_grids <- redundancy$id_grid[redundancy$redundancy == 0]
lych_cr <- lych_cr[!lych_cr$id_grid %in% inv_grids, ]

#=====================================================================================================#

#==============#
# PRUNING TREE #
#==============#

lych_pruned.tree <- drop.tip(lych_tree, 
                                   lych_tree$tip.label[!lych_tree$tip.label %in% lych_cr$gen_sp])

#==========================#
# RESCALING BRANCH LENGTHS #
#==========================#

edge.l <- c()
for(i in 1:length(lych_pruned.tree$edge.length)){
  edge.l[i] <- lych_pruned.tree$edge.length[i]/sum(lych_pruned.tree$edge.length)
}

lych_pruned.tree$edge.length <- edge.l

#=======#
# NEXUS #
#=======#

write.nexus(lych_pruned.tree, file = "trees/Lychnophorinae/pruned_tree-lych.nex")

#=====================================================================================================#

#========#
# MATRIX #
#========#

#Composing a presence/absence matrix
lych_matrix <- matrix(data = NA, nrow = length(unique(lych_cr$id_grid)), 
                            ncol = length(unique(lych_cr$gen_sp)))
lych_matrix <- as.data.frame(lych_matrix)
colnames(lych_matrix) <- unique(lych_cr$gen_sp)
rownames(lych_matrix) <- unique(lych_cr$id_grid)
for(i in 1:nrow(lych_matrix)){
  for(j in 1:ncol(lych_matrix)){
    if(colnames(lych_matrix)[j] %in% lych_cr$gen_sp[lych_cr$id_grid == rownames(lych_matrix)[i]]){
      lych_matrix[i, j] <- 1
    } else {
      lych_matrix[i, j] <- 0  
    }
  }
}

#Saving matrix
write.csv(lych_matrix, "matrices/lych_matrix.csv", row.names = TRUE)
