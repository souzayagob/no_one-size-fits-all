#=====================================================================================================#

library(rgdal)
library(recluster)
library(dendextend)
library(ape)
library(picante)

#======================================================================================================#

#=======#
# INPUT #
#=======#

#==========#
# Matrix #
#==========#

#Loading matrix
vel_matrix <- read.csv(file = "matrices/vel_matrix.csv", row.names = 1)

#Removing taxa that only have one recorded presence
vel_matrix <- vel_matrix[ , which(colSums(vel_matrix) > 1)]

#Removing empty sites 
vel_matrix <- vel_matrix[which(rowSums(vel_matrix) > 0), ]

#======#
# Tree #
#======#

#Loading tree
vel_tree <- read.nexus("trees/Velloziaceae/pruned_tree-vel.nex")

#Dropping tips (taxa that only have one recorded presence)
vel_tree <- drop.tip(vel_tree, 
                            vel_tree$tip.label[!vel_tree$tip.label %in% colnames(vel_matrix)])

#============#
# shapefiles #
#============#

#Loading cr grids and the Brazilian terrestrial territory
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp")
br <- readOGR("shapefiles/br_unidades_da_federacao/BRUFE250GC_SIR.shp") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#Projecting br
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
br <- spTransform(br, crswgs84)

#=======#
# Grids #
#=======#

#Data frame with all grid ids
grids_df <- as.data.frame(grids_cr@data)

#======================================================================================================#

#=========#
# UNIFRAC #
#=========#

#UniFrac matrix
unifrac_matrix <- as.matrix(unifrac(vel_matrix, vel_tree))
write.csv(file = "results/Velloziaceae/unifrac_matrix.csv", unifrac_matrix)

#===============#
# Running UPGMA #
#===============#

upgma <- recluster.cons(vel_matrix, vel_tree, dist = "unifrac",
                        tr = 1000, p = 0.5, method = "average")
upgma_cons <- upgma$cons
upgma_cons <- di2multi(upgma_cons) #identifying polytomies
hc <- as.hclust(upgma$cons) #dendrogram

#Fusion levels (useful to define the number of clusters)
plot(
  hc$height,
  nrow(vel_matrix):2,
  type = "S",
  main = "Fusion levels - Chord - UPGMA",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(hc$height,
     nrow(vel_matrix):2,
     nrow(vel_matrix):2,
     col = "red",
     cex = 0.8)

#====================#
# Dendogram and plot #
#====================#

#Number of clusters
ncluster <- 6
dend <- as.dendrogram(hc) 

#Defining colors
colors <- c("#41D91E",
            "#FB9A99",
            "#E300F7",
            "#FF0005",
            "#FF7F00",
            "#1F78B4",
            "#0000A3")

#Coloring clusters in the dendogram
dend <- color_branches(dend, k = ncluster, col = colors[1:ncluster])

#Indentifying the clusters of each grid and assigning them the respective color
groups_id <- data.frame("id" = labels(dend), "color" = get_leaves_branches_col(dend), 
                        "cluster_membership" = NA)
rb <- colors[1:ncluster]
names(rb) <- as.character(1:length(rb))
for(i in 1:nrow(groups_id)){
  groups_id$cluster_membership[i] <- names(rb)[unname(rb) == groups_id$color[i]]
}

#Merging everything into a polygon dataframe
unifrac_poly <- sp::merge(grids_cr, groups_id, by.x = "id")
unifrac_poly$cluster_membership <- factor(unifrac_poly$cluster_membership, 
                                          levels = unique(groups_id$cluster_membership))

#=====================================================================================================#

#=========#
# FIGURES #
#=========#

#Plot
unifrac_plot <- spplot(unifrac_poly, 
                       zcol = "cluster_membership", 
                       xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5), 
                       colorkey = TRUE, 
                       sp.layout = list(list(br, fill = "gray")), 
                       col.regions = colors[1:length(levels(unifrac_poly$cluster_membership))], 
                       scales = list(draw = FALSE))

#Dendogram
labels(dend) <- NULL
dend <- assign_values_to_branches_edgePar(dend = dend, value = 1, edgePar = "lwd")

pdf("plots/Velloziaceae/upgma/unifrac_plot-6.pdf",
    height = 4, width = 4); unifrac_plot; dev.off()
pdf("plots/Velloziaceae/upgma/unifrac_dend-6.pdf", 
    width = 4, height = 4); plot_horiz.dendrogram(dend, axes = F, center = TRUE); dev.off()