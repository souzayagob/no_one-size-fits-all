library(picante)
library(viridis)
library(tidyverse)
library(rgdal)
library(ape)
library(PDcalc) #https://rdrr.io/github/davidnipperess/PDcalc/man/phyloendemism.html

#=====================================================================================================#

#=======#
# INPUT #
#=======#

#Reading matrix
coman_matrix <- read.csv(file = "matrices/coman_matrix.csv", row.names = 1)

#Reading tree
coman_tree <- read.nexus("trees/Comanthera/pruned_tree-coman.nex")

#Reading grids
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Reading shapefile: Brazilian terrestrial territory
br <- readOGR("shapefiles/br_unidades_da_federacao/BRUFE250GC_SIR.shp") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#=====================================================================================================#

#====#
# PE #
#====#

#Calculating PE
pe <- data.frame("id" = rownames(coman_matrix), 
                 "pe" = phyloendemism(coman_matrix, coman_tree, weighted = T)) 

#PE randomization (independent swap)
rand.matrix.pe <- function(t, mat){
  phyloendemism(randomizeMatrix(mat, null.model = "independentswap"), t)[ , 1]
}
set.seed(7)
null.output <- replicate(999, rand.matrix.pe(coman_tree, coman_matrix))
rownames(null.output) <- rownames(coman_matrix)
ses.all <- (pe$pe - apply(null.output, MARGIN = 1, mean)) / apply(null.output,
                                                                  MARGIN = 1, sd)
p.val.all <- apply(cbind(pe$pe, null.output), MARGIN = 1, rank)[1, ] / 1000

#Concatenating results
pe$p.pe <- p.val.all
pe$ses.pe <- ses.all

#Significance
pe_poly <- merge(grids_cr, pe, by.x = "id")
pe_poly$significance <- NA
for(i in 1:nrow(pe_poly@data)){
  if(is.na(pe_poly$p.pe[i])){
    pe_poly$significance[i] <- NA
  } else if(pe_poly$p.pe[i] >= 0.99 ){
    pe_poly$significance[i] <-  "≥ 0.99"
  } else if(pe_poly$p.pe[i] >= 0.975 & pe_poly$p.pe[i] < 0.99){
    pe_poly$significance[i] <-  "≥ 0.975"
  } else if(pe_poly$p.pe[i] <= 0.025 & pe_poly$p.pe[i] > 0.01){
    pe_poly$significance[i] <-  "≥ 0.025"
  } else if(pe_poly$p.pe[i] <= 0.01){
    pe_poly$significance[i] <-  "≥ 0.01"
  } else {
    pe_poly$significance[i] <- "Not significant"
  }
}
pe_poly$significance <- factor(pe_poly$significance, levels = c("≥ 0.01", 
                                                                "≥ 0.025",
                                                                "Not significant",
                                                                "≥ 0.975",
                                                                "≥ 0.99"))

#=====================================================================================================#

#=========================#
# PE: tree for comparison #
#=========================#

#Tree for comparison
comp_tree <- coman_tree
comp_tree$edge.length <- rep(1/length(comp_tree$edge.length), 
                             length(comp_tree$edge.length)) #equal branch lengths

#Calculating PE
pe_comp <- data.frame("id" = rownames(coman_matrix), 
                      "pe_comp" = phyloendemism(coman_matrix, comp_tree, weighted = T)) 

#PE randomization (independent swap)
rand.matrix.pe <- function(t, mat){
  phyloendemism(randomizeMatrix(mat, null.model = "independentswap"), t)[ , 1]
}
set.seed(7)
null.output <- replicate(999, rand.matrix.pe(comp_tree, coman_matrix))
rownames(null.output) <- rownames(coman_matrix)
ses.all <- (pe_comp$pe_comp - apply(null.output, MARGIN=1, mean)) / apply(null.output,
                                                                          MARGIN = 1, sd)
p.val.all <- apply(cbind(pe_comp$pe_comp, null.output), MARGIN = 1, rank)[1, ] / 1000

#Concatenating results
pe_comp$p.pe_comp <- p.val.all
pe_comp$ses.pe_comp <- ses.all

#=====================================================================================================#

#=====#
# RPE #
#=====#

#Calculating RPE
rpe <-merge(pe, pe_comp)
rpe$rpe <- rpe$pe/rpe$pe_comp

#RPE randomization (independent swap)
rand.matrix.rpe <- function(t, mat){
  comp_t <- t
  comp_t$edge.length <- rep(1/length(comp_t$edge.length), 
                            length(comp_t$edge.length))
  rand.mat <- randomizeMatrix(mat, null.model = "independentswap")
  phyloe <- phyloendemism(rand.mat, t)[ , 1]
  phyloe_comp <- phyloendemism(rand.mat, comp_t)[ , 1]
  phyloe/phyloe_comp  
}
set.seed(7)
null.output <- replicate(999, rand.matrix.rpe(coman_tree, coman_matrix))
rownames(null.output) <- rownames(coman_matrix)
ses.all <- (rpe$rpe - apply(null.output, MARGIN=1, mean)) / apply(null.output,
                                                                  MARGIN = 1, sd)
p.val.all <- apply(cbind(rpe$rpe, null.output), MARGIN = 1, rank)[1, ] / 1000

#Concatenating results
rpe$p.rpe <- p.val.all
rpe$ses.rpe <- ses.all

#Significance
rpe_poly <- merge(grids_cr, rpe, by.x = "id")
rpe_poly$significance <- NA
for(i in 1:nrow(rpe_poly@data)){
  if(is.na(rpe_poly$p.rpe[i])){
    rpe_poly$significance[i] <- NA
  } else if(rpe_poly$p.rpe[i] >= 0.99 ){
    rpe_poly$significance[i] <-  "≥ 0.99"
  } else if(rpe_poly$p.rpe[i] >= 0.975 & rpe_poly$p.rpe[i] < 0.99){
    rpe_poly$significance[i] <-  "≥ 0.975"
  } else if(rpe_poly$p.rpe[i] <= 0.025 & rpe_poly$p.rpe[i] > 0.01){
    rpe_poly$significance[i] <-  "≥ 0.025"
  } else if(rpe_poly$p.rpe[i] <= 0.01){
    rpe_poly$significance[i] <-  "≥ 0.01"
  } else {
    rpe_poly$significance[i] <- "Not significant"
  }
}
rpe_poly$significance <- factor(rpe_poly$significance, levels = c("≥ 0.01", 
                                                                  "≥ 0.025",
                                                                  "Not significant",
                                                                  "≥ 0.975",
                                                                  "≥ 0.99"))

#=====================================================================================================#

#========#
# CANAPE #
#========#

canape <- rpe
canape$step_one <- c()
for(i in 1:nrow(canape)){
  if(canape$p.pe[i] >= 0.95 | canape$p.pe_comp[i] >= 0.95){
    canape$step_one[i] <- TRUE
  } else{
    canape$step_one[i] <- FALSE
  }
}

canape$results <- NULL
for(i in 1:nrow(canape)){
  if(canape$step_one[i] == FALSE){
    canape$results[i] <- "Not significant"
  } else if(canape$p.rpe[i] >= 0.975 & canape$step_one[i] == TRUE){
    canape$results[i] <- "Paleo-endemism"
  } else if(canape$p.rpe[i] <= 0.025 & canape$step_one[i] == TRUE){
    canape$results[i] <- "Neo-endemism"
  } else if(canape$p.pe[i] >= 0.95 & canape$p.pe_comp[i] >= 0.95 & canape$step_one[i] == TRUE){
    canape$results[i] <- "Mixed endemism"
  }else{
    canape$results[i] <- "Uncertain" #This should be checked later. 
  }
  if(canape$p.pe[i] >= 0.99 & canape$p.pe_comp[i] >= 0.99 & canape$step_one[i] == TRUE){
    canape$results[i] <- "Super endemism"
  }
} 

#Plotting
canape_poly <- merge(grids_cr, canape[ , which(colnames(canape) %in% c("id", "results"))], by.x = "id")
canape_poly$results <- factor(canape_poly$results, levels = c("Neo-endemism",
                                                              "Paleo-endemism",
                                                              "Not significant",
                                                              "Uncertain",
                                                              "Mixed endemism",
                                                              "Super endemism"))

#=====================================================================================================#

#=========#
# FIGURES #
#=========#

#PE
pe_plot <- spplot(pe_poly,
                  zcol = "pe",
                  xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                  par.settings=list(fontsize = list(text = 18)),
                  colorkey = TRUE, 
                  sp.layout = list(list(br, fill = "gray")), 
                  col.regions = rev(magma(16)), 
                  scales = list(draw = FALSE))

#PE significance
peP_plot <- spplot(pe_poly,  
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

#RPE
rpe_plot <- spplot(rpe_poly,
                   zcol = "rpe",
                   xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                   par.settings=list(fontsize = list(text = 18)),
                   colorkey = TRUE, 
                   sp.layout = list(list(br, fill = "gray")), 
                   col.regions = rev(magma(16)), 
                   scales = list(draw = FALSE))

#RPE #significance
rpeP_plot <- spplot(rpe_poly,  
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

#CANAPE
canape_plot <- spplot(canape_poly, zcol = "results", colorkey = TRUE, 
                      sp.layout = list(list(br, fill = "gray")), 
                      xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                      col.regions = c("#231151FF",
                                      "#B63679FF",
                                      "khaki1",
                                      "#FFFFCC",
                                      "#4AC16DFF",
                                      "#CFE11CFF"), scales = list(draw = FALSE))

png("plots/Comanthera/canape/pe_plot.png",
    height = 4, width = 4, units = 'in', res=300); pe_plot; dev.off()
png("plots/Comanthera/canape/peP_plot.png",
    height = 4, width = 4, units = 'in', res=300); peP_plot; dev.off()
png("plots/Comanthera/canape/rpe_plot.png",
    height = 4, width = 4, units = 'in', res=300); rpe_plot; dev.off()
png("plots/Comanthera/canape/rpeP_plot.png",
    height = 4, width = 4, units = 'in', res=300); rpeP_plot; dev.off()
png("plots/Comanthera/canape/canape_plot.png",
    height = 4, width = 4, units = 'in', res=300); canape_plot; dev.off()
