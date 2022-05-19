#======================================================================================================#

library(rgdal)
library(recluster)
library(vegan)
#library(phytools)
#library(maps)
#library(stats)
#library(cluster)
library(dendextend)

#======================================================================================================#

#=======#
# INPUT #
#=======#

#==========#
# Matrices #
#==========#

#Mimosa
mimosa_matrix <- read.csv(file = "matrices/mimosa_matrix.csv", row.names = 1)

#Diplusodon
diplusodon_matrix <- read.csv(file = "matrices/diplusodon_matrix.csv", row.names = 1)

#Comanthera
coman_matrix <- read.csv(file = "matrices/coman_matrix.csv", row.names = 1)

#Velloziaceae
vel_matrix <- read.csv(file = "matrices/vel_matrix.csv", row.names = 1)

#Lychnophorinae
lych_matrix <- read.csv(file = "matrices/lych_matrix.csv", row.names = 1)

#Paepalanthus
paep_matrix <- read.csv(file = "matrices/paep_matrix.csv", row.names = 1)

#Calliandra
calli_matrix <- read.csv(file = "matrices/paep_matrix.csv", row.names = 1)

#Minaria
minaria_matrix <- read.csv(file = "matrices/paep_matrix.csv", row.names = 1)

#Marcetieae
marcetieae_matrix <- read.csv(file = "matrices/paep_matrix.csv", row.names = 1)

#============#
# shapefiles #
#============#

#Loading cr grids and the Brazilian terrestrial territory
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp")
br <- readOGR("shapefiles/br_unidades_da_federacao/BRUFE250GC_SIR.shp") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#Projecting br
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
br <- spTransform(br, crswgs84)

#=============#
# Data frames #
#=============#

#Data frame with all grid ids
grids_df <- as.data.frame(grids_cr@data)

#Matrix - synthesis
matrices <- ls()[grep("matrix", ls())]
for(i in 1:length(matrices)){
  m <- get(matrices[i])
  m$id <- rownames(m)
  assign(matrices[i], m)
  rm(m)
}
all_m <- plyr::join(grids_df, mimosa_matrix, by = "id")
all_m <- plyr::join(all_m[ , -which(colnames(all_m) == "layer")], diplusodon_matrix, by = "id")
all_m <- plyr::join(all_m, coman_matrix, by = "id")
all_m <- plyr::join(all_m, vel_matrix, by = "id")
all_m <- plyr::join(all_m, lych_matrix, by = "id")
all_m <- plyr::join(all_m, paep_matrix, by = "id")
all_m <- plyr::join(all_m, calli_matrix, by = "id")
all_m <- plyr::join(all_m, minaria_matrix, by = "id")
all_m <- plyr::join(all_m, marcetieae_matrix, by = "id")

#Replacing NAs by 0
all_m[is.na(all_m)] <- 0

#Renaming rows
rownames(all_m) <- all_m$id

#Removing id column
all_m <- all_m[ , -which(colnames(all_m) == "id")]

#Removing taxa that only have one recorded presence
all_m <- all_m[ , which(colSums(all_m) > 1)]

#Removing empty sites 
all_m <- all_m[which(rowSums(all_m) > 0), ]

#======================================================================================================#

#=========#
# JACCARD #
#=========#

#Jaccard distance matrix
jaccard_matrix <- as.matrix(vegdist(all_m, method = "jaccard", diag = TRUE))
write.csv(file = "results/synthesis/jaccard_matrix.csv", jaccard_matrix)

#===============#
# Running UPGMA #
#===============#

upgma <- recluster.cons(all_m, dist = "jaccard",
                        tr = 1000, p = 0.5, method = "average")
upgma_cons <- upgma$cons
upgma_cons <- di2multi(upgma_cons) #identifying polytomies
hc <- as.hclust(upgma$cons) #dendrogram

#Fusion levels (useful to define the number of clusters)
plot(
  hc$height,
  nrow(all_m):2,
  type = "S",
  main = "Fusion levels - Chord - UPGMA",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(hc$height,
     nrow(all_m):2,
     nrow(all_m):2,
     col = "red",
     cex = 0.8)

#====================#
# Dendogram and plot #
#====================#

#Number of clusters
ncluster <- 14
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
                       colorkey = TRUE, sp.layout = list(list(br, fill = "grey")), 
                       col.regions = colors[1:length(levels(jaccard_poly$cluster_membership))], 
                       scales = list(draw = FALSE))

#Dendogram
labels(dend) <- NULL
dend <- assign_values_to_branches_edgePar(dend = dend, value = 4, edgePar = "lwd")

png("plots/synthesis/jaccard/jaccard_plot.png",
    height = 4, width = 4, units = 'in', res=300); jaccard_plot; dev.off()
cairo_pdf("plots/synthesis/jaccard/jaccard_dend.pdf", 
          width = 11, height = 11); plot_horiz.dendrogram(dend, axes = F); dev.off()
