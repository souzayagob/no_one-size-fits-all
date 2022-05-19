#=====================================================================================================#

library(picante)
library(viridis)
library(colorspace)
library(tidyverse)
library(rgdal)
library(ape)

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading matrix
vel_matrix <- read.csv(file = "matrices/vel_matrix.csv", row.names = 1)

#Reading tree
vel_tree <- read.nexus("trees/Velloziaceae/pruned_tree-vel.nex")

#Reading grids
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Reading shapefile: Brazilian terrestrial territory
br <- readOGR("shapefiles/br_unidades_da_federacao/BRUFE250GC_SIR.shp") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#=====================================================================================================#

#====#
# PD #
#====#

#Setting random number generator and running phylogenetic diversity analysis
set.seed(7)
pd_stats <- ses.pd(vel_matrix, vel_tree, include.root = TRUE, null.model = "taxa.label")

#Creating a column specifying grid ids
pd_stats$id <- rownames(pd_stats)

#Running a linear regression and incorporating residuals
pd_stats$residuals <- lm(pd.obs ~ ntaxa, pd_stats)$res

#Rounding values in order to standardize plot size (except for p-values, id and ntaxa)
pd_stats[ , 
          -which(colnames(pd_stats) %in% c("ntaxa",
                                           "pd.obs.p", 
                                           "id"))] <- round(pd_stats[ , 
                                                                      -which(colnames(pd_stats) %in% c("ntaxa",
                                                                                                       "pd.obs.p", 
                                                                                                       "id"))], 
                                                            2)

#Saving results
#write.csv(file = "results/Velloziaceae/pd_stats.csv", pd_stats, row.names = FALSE)

#Loading results
#pd_stats <- read.csv("results/Velloziaceae/pd_stats.csv")

#Assigning values for each grid by merging results and spatial grids
pd_poly <- merge(grids_cr, pd_stats, by.x = "id")

#Defining cells with statistical significant PD
pd_poly$significance <- NA
for(i in 1:nrow(pd_poly@data)){
  if(is.na(pd_poly$pd.obs.p[i])){
    pd_poly$significance[i] <- NA
  } else if(pd_poly$pd.obs.p[i] >= 0.99 ){
    pd_poly$significance[i] <-  "≥ 0.99"
  } else if(pd_poly$pd.obs.p[i] >= 0.975 & pd_poly$pd.obs.p[i] < 0.99){
    pd_poly$significance[i] <-  "≥ 0.975"
  } else if(pd_poly$pd.obs.p[i] <= 0.025 & pd_poly$pd.obs.p[i] > 0.01){
    pd_poly$significance[i] <-  "≥ 0.025"
  } else if(pd_poly$pd.obs.p[i] <= 0.01){
    pd_poly$significance[i] <-  "≥ 0.01"
  } else {
    pd_poly$significance[i] <- "Not significant"
  }
}
pd_poly$significance <- factor(pd_poly$significance, levels = c("≥ 0.01", 
                                                                "≥ 0.025",
                                                                "Not significant",
                                                                "≥ 0.975",
                                                                "≥ 0.99"))

#=====================================================================================================#

#=========#
# FIGURES #
#=========#

#Labels limits based on the range of each attribute. This was retrieved using summary(pd_poly@data)
summary(pd_poly@data)

#PD
pd_plot <- spplot(pd_poly,
                  zcol = "pd.obs",
                  xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                  par.settings=list(fontsize = list(text = 21)),
                  colorkey = TRUE, 
                  at = seq(0, 0.561, length.out = 16),
                  sp.layout = list(list(br, fill = "gray")), 
                  col.regions = rev(magma(16)), 
                  scales = list(draw = FALSE))

#SR
sr_plot <- spplot(pd_poly,
                  zcol = "ntaxa", 
                  xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                  par.settings=list(fontsize = list(text = 21)),
                  at = seq(0, 51, length.out = 16), 
                  sp.layout = list(list(br, fill = "gray")), 
                  col.regions = rev(magma(16)), 
                  scales = list(draw = FALSE))

#Correlation: PD and SR
corPdsr <- as.character(formatC(cor(pd_stats$ntaxa, 
                                    pd_stats$pd.obs, use = "na.or.complete"))) 
corPdsr_plot <- ggplot(data = pd_stats, mapping = aes(jitter(ntaxa), pd.obs))+
  geom_jitter()+
  labs(title = "Velloziaceae", subtitle = paste("r =", corPdsr))+
  xlab("Richness")+
  ylab("Phylogenetic diversity")+
  theme_bw(base_size = 14)

#SES
pdses_plot <- spplot(pd_poly,
                     zcol = "pd.obs.z",
                     xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                     par.settings=list(fontsize = list(text = 21)),
                     at = seq(-4.451, 4.451, length.out = 16),
                     colorkey = TRUE, 
                     sp.layout = list(list(br, fill = "gray")), 
                     col.regions = viridis(16), 
                     scales = list(draw = FALSE))

#Residuals
pdres_plot <- spplot(pd_poly,
                     zcol = "residuals",
                     xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                     par.settings=list(fontsize = list(text = 21)),
                     at = seq(-0.101, 0.101, length.out = 16),
                     colorkey = TRUE, 
                     sp.layout = list(list(br, fill = "gray")), 
                     col.regions = c(viridis(16)[1:7], "khaki1", rev(heat.colors(16)[1:8])), 
                     scales = list(draw = FALSE))

#PD statistical significance
pdP_plot <- spplot(pd_poly,  
                   zcol = "significance", 
                   xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                   colorkey = TRUE, 
                   sp.layout = list(list(br, fill = "gray")), 
                   col.regions = c("firebrick4",
                                   "firebrick3",
                                   "khaki1",
                                   "purple3",
                                   "purple4"), 
                   scales = list(draw = FALSE))

#Saving plots
png("plots/Velloziaceae/pd/pd_plot.png",
    height = 4, width = 4, units = 'in', res=300); pd_plot; dev.off()
png("plots/Velloziaceae/pd/sr_plot.png",
    height = 4, width = 4, units = 'in', res=300); sr_plot; dev.off()
png("plots/Velloziaceae/pd/corPdsr_plot.png",
    height = 4, width = 4, units = 'in', res=300); corPdsr_plot; dev.off()
png("plots/Velloziaceae/pd/pdses_plot.png",
    height = 4, width = 4, units = 'in', res=300); pdses_plot; dev.off()
png("plots/Velloziaceae/pd/pdres_plot.png",
    height = 4, width = 4, units = 'in', res=300); pdres_plot; dev.off()
png("plots/Velloziaceae/pd/pdP_plot.png",
    height = 4, width = 4, units = 'in', res=300); pdP_plot; dev.off()
