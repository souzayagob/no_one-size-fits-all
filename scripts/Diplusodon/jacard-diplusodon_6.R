#=====================================================================================================#

library(rgdal)
library(recluster)
library(vegan)
library(dendextend)

#======================================================================================================#

#=======#
# INPUT #
#=======#

#==========#
# Matrix #
#==========#

#Loading matrix
diplusodon_matrix <- read.csv(file = "matrices/diplusodon_matrix.csv", row.names = 1)

#Removing taxa that only have one recorded presence
diplusodon_matrix <- diplusodon_matrix[ , which(colSums(diplusodon_matrix) > 1)]

#Removing empty sites 
diplusodon_matrix <- diplusodon_matrix[which(rowSums(diplusodon_matrix) > 0), ]

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
# JACCARD #
#=========#

#Jaccard distance matrix
jaccard_matrix <- as.matrix(vegdist(diplusodon_matrix, method = "jaccard", diag = TRUE))
write.csv(file = "results/Diplusodon/jaccard_matrix.csv", jaccard_matrix)

#CODE FOR CORRELATION
#jaccard_matrix <- as.matrix(read.csv("results/Diplusodon/jaccard_matrix.csv", row.names = 1))
#colnames(jaccard_matrix) <- gsub("X", "", colnames(jaccard_matrix))
#unfirac_matrix <- as.matrix(read.csv("results/Diplusodon/unifrac_matrix.csv", row.names = 1))
#colnames(unifrac_matrix) <- gsub("X", "", colnames(unifrac_matrix))
#cor(c(jaccard_matrix), c(unifrac_matrix))

#===============#
# Running UPGMA #
#===============#

upgma <- recluster.cons(diplusodon_matrix, dist = "jaccard",
                        tr = 1000, p = 0.5, method = "average")
upgma_cons <- upgma$cons
upgma_cons <- di2multi(upgma_cons) #identifying polytomies
hc <- as.hclust(upgma$cons) #dendrogram

#Fusion levels (useful to define the number of clusters)
plot(
  hc$height,
  nrow(diplusodon_matrix):2,
  type = "S",
  main = "Fusion levels - Chord - UPGMA",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(hc$height,
     nrow(diplusodon_matrix):2,
     nrow(diplusodon_matrix):2,
     col = "red",
     cex = 0.8)

#====================#
# Dendrogram and plot #
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
dend <- color_branches(dend, k = ncluster, col = colors[1:ncluster], groupLabels = F)

#Indentifying the clusters of each grid and assigning them the respective color
groups_id <- data.frame("id" = labels(dend), "color" = get_leaves_branches_col(dend), 
                        "cluster_membership" = NA)
rb <- colors[1:ncluster]
names(rb) <- as.character(1:length(rb))
for(i in 1:nrow(groups_id)){
  groups_id$cluster_membership[i] <- names(rb)[unname(rb) == groups_id$color[i]]
}

#Merging everything into a polygon dataframe
jaccard_poly <- sp::merge(grids_cr, groups_id, by.x = "id")
jaccard_poly$cluster_membership <- factor(jaccard_poly$cluster_membership, 
                                          levels = unique(groups_id$cluster_membership))

#=====================================================================================================#

#=========#
# FIGURES #
#=========#

#Plot
jaccard_plot <- spplot(jaccard_poly, 
                       zcol = "cluster_membership", 
                       xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5), 
                       colorkey = FALSE, 
                       sp.layout = list(list(br, fill = "grey90")), 
                       col.regions = colors[1:length(levels(jaccard_poly$cluster_membership))], 
                       scales = list(draw = FALSE))

#Dendrogram
labels(dend) <- NULL
dend <- assign_values_to_branches_edgePar(dend = dend, value = 1, edgePar = "lwd")

pdf("plots/Diplusodon/upgma/jaccard_plot-6.pdf",
    height = 4, width = 4); jaccard_plot; dev.off()
pdf("plots/Diplusodon/upgma/jaccard_dend-6.pdf", 
    width = 4, height = 4); plot_horiz.dendrogram(dend, axes = F, center = TRUE); dev.off()