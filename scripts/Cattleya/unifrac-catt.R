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
catt_matrix <- read.csv(file = "matrices/catt_matrix.csv", row.names = 1)

#Removing taxa that only have one recorded presence
catt_matrix <- catt_matrix[ , which(colSums(catt_matrix) > 1)]

#Removing empty sites 
catt_matrix <- catt_matrix[which(rowSums(catt_matrix) > 0), ]

#======#
# Tree #
#======#

#Loading tree
catt_tree <- read.nexus("trees/Cattleya/pruned_tree-catt.nex")

#Dropping tips (taxa that only have one recorded presence)
catt_tree <- drop.tip(catt_tree, 
                       catt_tree$tip.label[!catt_tree$tip.label %in% colnames(catt_matrix)])

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
unifrac_matrix <- as.matrix(unifrac(catt_matrix, catt_tree))
write.csv(file = "results/Cattleya/unifrac_matrix.csv", unifrac_matrix)

#===============#
# Running UPGMA #
#===============#

upgma <- recluster.cons(catt_matrix, catt_tree, dist = "unifrac",
                        tr = 1000, p = 0.5, method = "average")
upgma_cons <- upgma$cons
upgma_cons <- di2multi(upgma_cons) #identifying polytomies
hc <- as.hclust(upgma$cons) #dendrogram

#Fusion levels (useful to define the number of clusters)
plot(
  hc$height,
  nrow(catt_matrix):2,
  type = "S",
  main = "Fusion levels - Chord - UPGMA",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(hc$height,
     nrow(catt_matrix):2,
     nrow(catt_matrix):2,
     col = "red",
     cex = 0.8)

#====================#
# Dendogram and plot #
#====================#

#Number of clusters
ncluster <- 10
dend <- as.dendrogram(hc) 

#Defining colors
colors <- c("#41D91E", 
            "#FB9A99",
            "#003200",
            "#008805",
            "#E300F7",
            "#FFB559",
            "#FF0005",
            "#FF7F00",
            "#6A3D9A",
            "#1F78B4",
            "#FFFF99",
            "#0000A3",
            "#FAD900",
            "#B15928")

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
dend <- assign_values_to_branches_edgePar(dend = dend, value = 4, edgePar = "lwd")

png("plots/Cattleya/upgma/unifrac_plot.png",
    height = 4, width = 4, units = 'in', res=300); unifrac_plot; dev.off()
cairo_pdf("plots/Cattleya/upgma/unifrac_dend.pdf", 
          width = 11, height = 11); plot_horiz.dendrogram(dend, axes = F); dev.off()