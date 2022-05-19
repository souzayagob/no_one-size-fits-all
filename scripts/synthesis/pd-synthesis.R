#=====================================================================================================#

library(picante)
library(viridis)
library(tidyverse)
library(rgdal)
library(ape)

#=====================================================================================================#

#=======#    
# INPUT #
#=======#

#========#
# Mimosa #
#========#

#Loading PD results
mimosa_pd <- read.csv(file = "results/Mimosa/pd_stats.csv")

#============#
# Diplusodon #
#============#

#Loading PD results
diplusodon_pd <- read.csv(file = "results/Diplusodon/pd_stats.csv")

#============#
# Comanthera #
#============#

#Loading PD results
coman_pd <- read.csv(file = "results/Comanthera/pd_stats.csv")

#==============#
# Velloziaceae #
#==============#

#Loading PD results
vel_pd <- read.csv(file = "results/Velloziaceae/pd_stats.csv")

#================#
# Lychnophorinae #
#================#

#Loading PD results
lych_pd <- read.csv(file = "results/Lychnophorinae/pd_stats.csv")

#==============#
# Paepalanthus #
#==============#

#Loading PD results
paep_pd <- read.csv(file = "results/Paepalanthus/pd_stats.csv")

#============#
# Calliandra #
#============#

#Loading PD results
calli_pd <- read.csv(file = "results/Calliandra/pd_stats.csv")

#=========#
# Minaria #
#=========#

#Loading PD results
minaria_pd <- read.csv(file = "results/Minaria/pd_stats.csv")

#===========#
# Habenaria #
#===========#

#Loading PD results
#habenaria_pd <- read.csv(file = "results/Habenaria/pd_stats.csv")

#==========#
# Cattleya #
#==========#

#Loading PD results
#catt_pd <- read.csv(file = "results/Cattleya/pd_stats.csv")


#============#
# Trimezieae #
#============#

#Loading PD results
#trimezieae_pd <- read.csv(file = "results/Trimezieae/pd_stats.csv")

#============#
# Marcetieae #
#============#

#Loading PD results
marcetieae_pd <- read.csv(file = "results/Marcetieae/pd_stats.csv")


#=====================#
# Weighting SR and PD #
#=====================#

results <- ls()
for(i in 1:length(results)){
  df <- get(results[i])
  df$ntaxa_w <- df$ntaxa/sum(df$ntaxa)
  df$pd.obs_w <- df$pd.obs/sum(df$pd.obs)
  assign(results[i], df)
  rm(df)
}

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

#Data frame containing sum of SR and PD for each ID. Weighted (conc_w) and non-weigthted (conc)
conc <- plyr::join(grids_df, mimosa_pd, by = "id")
conc <- plyr::join(conc[ , -which(colnames(conc) == "layer")], diplusodon_pd, by = "id")
conc <- plyr::join(conc, coman_pd, by = "id")
conc <- plyr::join(conc, vel_pd, by = "id")
conc <- plyr::join(conc, lych_pd, by = "id")
conc <- plyr::join(conc, paep_pd, by = "id")
conc <- plyr::join(conc, calli_pd, by = "id")
conc <- plyr::join(conc, minaria_pd, by = "id")
#conc <- plyr::join(conc, habenaria_pd, by = "id")
#conc <- plyr::join(conc, catt_pd, by = "id")
#conc <- plyr::join(conc, trimezieae_pd, by = "id")
conc <- plyr::join(conc, marcetieae_pd, by = "id")
conc <- as.data.frame(sapply(split.default(conc, names(conc)), rowSums, na.rm = TRUE))
conc_w <- conc[-which(conc$ntaxa == 0), which(colnames(conc) %in% c("id", 
                                                                  "ntaxa_w",
                                                                  "pd.obs_w"))]
conc <- conc[-which(conc$ntaxa == 0), which(colnames(conc) %in% c("id", 
                                                                  "ntaxa",
                                                                  "pd.obs"))]

#=====================================================================================================#

#===========#
# SYNTHESIS #
#===========#

#Linear model and residuals
set.seed(7)
pd_lm.total <- lm(pd.obs ~ ntaxa, conc)
conc$total_residuals <- pd_lm.total$residuals

#Merging synthesis with spatial grids
pdPoly <- merge(grids_cr, conc, by.x = "id")

#Rounding values in order to standardize plot size (id and ntaxa)
pdPoly@data[ , 
          -which(colnames(pdPoly@data) %in% c("ntaxa",
                                           "id"))] <- round(pdPoly@data[ , 
                                                                      -which(colnames(pdPoly@data) %in% c("ntaxa",
                                                                                                       "id"))], 2)

#=======#
# Plots #
#=======#

#Labels limits based on the range of each attribute. This was retrieved using summary(pd_poly@data)
summary(pdPoly@data)

#SR
sr_plot <- spplot(pdPoly,
                   zcol = "ntaxa",
                   xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                   par.settings=list(fontsize = list(text = 21)),
                   at = seq(0, 224, length.out = 16),
                   colorkey = TRUE, sp.layout = list(list(br, 
                                                          fill = "gray")), 
                   col.regions = rev(magma(16)), scales = list(draw = FALSE))
#png("plots/synthesis/pd/sr_plot.png",
#    height = 4, width = 4, units = 'in', res=300); sr_plot; dev.off()

#PD
pd_plot <- spplot(pdPoly,
                   zcol = "pd.obs",
                   xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                   par.settings=list(fontsize = list(text = 21)),
                   at = seq(0, 4.851, length.out = 16),
                   colorkey = TRUE, sp.layout = list(list(br, 
                                                          fill = "gray")), 
                   col.regions = rev(magma(16)), scales = list(draw = FALSE))
#png("plots/synthesis/pd/pd_plot.png",
#    height = 4, width = 4, units = 'in', res=300); pd_plot; dev.off()

#Residuals
pdres_plot <- spplot(pdPoly,
                     zcol = "total_residuals",
                     xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                     at = seq(-1.031, 1.031, length.out = 16),
                     par.settings=list(fontsize = list(text = 21)),
                     colorkey = TRUE, sp.layout = list(list(br, 
                                                            fill = "gray")), 
                     col.regions = c(viridis(16)[1:7], "khaki1", rev(heat.colors(16)[1:8])), scales = list(draw = FALSE))
#png("plots/synthesis/pd/pdres_plot.png",
#    height = 4, width = 4, units = 'in', res=300); pdres_plot; dev.off()

#==============#
# CORRELATIONS #
#==============#

#Joining concatenated dataset and individual datasets
corr_mimosa <- plyr::join(conc[ , c("id", "total_residuals")], mimosa_pd, by = "id")
corr_diplusodon <- plyr::join(conc[ , c("id", "total_residuals")], diplusodon_pd, by = "id")
corr_coman <- plyr::join(conc[ , c("id", "total_residuals")], coman_pd, by = "id")
corr_vel <- plyr::join(conc[ , c("id", "total_residuals")], vel_pd, by = "id")
corr_lych <- plyr::join(conc[ , c("id", "total_residuals")], lych_pd, by = "id")
corr_paep <- plyr::join(conc[ , c("id", "total_residuals")], paep_pd, by = "id")
corr_calli <- plyr::join(conc[ , c("id", "total_residuals")], calli_pd, by = "id")
corr_minaria <- plyr::join(conc[ , c("id", "total_residuals")], minaria_pd, by = "id")
#corr_habenaria <- plyr::join(conc[ , c("id", "total_residuals")], habenaria_pd, by = "id")
#corr_catt <- plyr::join(conc[ , c("id", "total_residuals")], catt_pd, by = "id")
#corr_trimezieae <- plyr::join(conc[ , c("id", "total_residuals")], trimezieae_pd, by = "id")
corr_marcetieae <- plyr::join(conc[ , c("id", "total_residuals")], marcetieae_pd, by = "id")


#=======#
# Plots #
#=======#

#Mimosa
mimosa_corrRes <- as.character(formatC(cor(corr_mimosa$residuals, 
                                    corr_mimosa$total_residuals, use = "na.or.complete"))) 
mimosa_corrRes.plot <- ggplot(data = corr_mimosa, mapping = aes(residuals, total_residuals))+
  geom_point(position = "jitter")+
  labs(title = "Mimosa", subtitle = paste("r =", mimosa_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))
png("plots/synthesis/pd/res_corr/mimosa.png",
    height = 4, width = 4, units = 'in', res=300); mimosa_corrRes.plot; dev.off()

#Diplusodon
diplusodon_corrRes <- as.character(formatC(cor(corr_diplusodon$residuals, 
                                          corr_diplusodon$total_residuals, use = "na.or.complete"))) 
diplusodon_corrRes.plot <- ggplot(data = corr_diplusodon, mapping = aes(residuals, total_residuals))+
  geom_point(position = "jitter")+
  labs(title = "Diplusodon", subtitle = paste("r =", diplusodon_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))
png("plots/synthesis/pd/res_corr/diplusodon.png",
    height = 4, width = 4, units = 'in', res=300); diplusodon_corrRes.plot; dev.off()

#Comanthera
coman_corrRes <- as.character(formatC(cor(corr_coman$residuals, 
                                          corr_coman$total_residuals, use = "na.or.complete"))) 
coman_corrRes.plot <- ggplot(data = corr_coman, mapping = aes(residuals, total_residuals))+
  geom_point(position = "jitter")+
  labs(title = "Comanthera", subtitle = paste("r =", coman_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))
png("plots/synthesis/pd/res_corr/coman.png",
    height = 4, width = 4, units = 'in', res=300); coman_corrRes.plot; dev.off()

#Velloziaceae
vel_corrRes <- as.character(formatC(cor(corr_vel$residuals, 
                                          corr_vel$total_residuals, use = "na.or.complete"))) 
vel_corrRes.plot <- ggplot(data = corr_vel, mapping = aes(residuals, total_residuals))+
  geom_point(position = "jitter")+
  labs(title = "Velloziaceae", subtitle = paste("r =", vel_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw(base_size = 14)+
  theme()
png("plots/synthesis/pd/res_corr/vel.png",
    height = 4, width = 4, units = 'in', res=300); vel_corrRes.plot; dev.off()

#Lychnophorinae
lych_corrRes <- as.character(formatC(cor(corr_lych$residuals, 
                                          corr_lych$total_residuals, use = "na.or.complete"))) 
lych_corrRes.plot <- ggplot(data = corr_lych, mapping = aes(residuals, total_residuals))+
  geom_point(position = "jitter")+
  labs(title = "Lychnophorinae", subtitle = paste("r =", lych_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw(base_size = 14)+
  theme()
png("plots/synthesis/pd/res_corr/lych.png",
    height = 4, width = 4, units = 'in', res=300); lych_corrRes.plot; dev.off()

#Paepalanthus
paep_corrRes <- as.character(formatC(cor(corr_paep$residuals, 
                                          corr_paep$total_residuals, use = "na.or.complete"))) 
paep_corrRes.plot <- ggplot(data = corr_paep, mapping = aes(residuals, total_residuals))+
  geom_point(position = "jitter")+
  labs(title = "Paepalanthus", subtitle = paste("r =", paep_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))
png("plots/synthesis/pd/res_corr/paep.png",
    height = 4, width = 4, units = 'in', res=300); paep_corrRes.plot; dev.off()

#Calliandra
calli_corrRes <- as.character(formatC(cor(corr_calli$residuals, 
                                         corr_calli$total_residuals, use = "na.or.complete"))) 
calli_corrRes.plot <- ggplot(data = corr_calli, mapping = aes(residuals, total_residuals))+
  geom_point(position = "jitter")+
  labs(title = "Calliandra", subtitle = paste("r =", calli_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))
png("plots/synthesis/pd/res_corr/calli.png",
    height = 4, width = 4, units = 'in', res=300); calli_corrRes.plot; dev.off()

#Minaria
minaria_corrRes <- as.character(formatC(cor(corr_minaria$residuals, 
                                         corr_minaria$total_residuals, use = "na.or.complete"))) 
minaria_corrRes.plot <- ggplot(data = corr_minaria, mapping = aes(residuals, total_residuals))+
  geom_point(position = "jitter")+
  labs(title = "Minaria", subtitle = paste("r =", minaria_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))
png("plots/synthesis/pd/res_corr/minaria.png",
    height = 4, width = 4, units = 'in', res=300); minaria_corrRes.plot; dev.off()

#Habenaria
#habenaria_corrRes <- as.character(formatC(cor(corr_habenaria$residuals, 
#                                         corr_habenaria$total_residuals, use = "na.or.complete"))) 
#habenaria_corrRes.plot <- ggplot(data = corr_habenaria, mapping = aes(residuals, total_residuals))+
#  geom_point(position = "jitter")+
#  labs(title = "Habenaria", subtitle = paste("r =", habenaria_corrRes))+
#  xlab("Residuals")+
#  ylab("Total residuals")+
#  theme_bw(base_size = 14)+
#  theme(plot.title = element_text(face = "italic"))
#png("plots/synthesis/pd/res_corr/habenaria.png",
#    height = 4, width = 4, units = 'in', res=300); habenaria_corrRes.plot; dev.off()

#Cattleya
#catt_corrRes <- as.character(formatC(cor(corr_catt$residuals, 
#                                         corr_catt$total_residuals, use = "na.or.complete"))) 
#catt_corrRes.plot <- ggplot(data = corr_catt, mapping = aes(residuals, total_residuals))+
#  geom_point(position = "jitter")+
#  labs(title = "Cattleya", subtitle = paste("r =", catt_corrRes))+
#  xlab("Residuals")+
#  ylab("Total residuals")+
#  theme_bw(base_size = 14)+
#  theme(plot.title = element_text(face = "italic"))
#png("plots/synthesis/pd/res_corr/catt.png",
#    height = 4, width = 4, units = 'in', res=300); catt_corrRes.plot; dev.off()

#Trimezieae
#trimezieae_corrRes <- as.character(formatC(cor(corr_trimezieae$residuals, 
#                                         corr_trimezieae$total_residuals, use = "na.or.complete"))) 
#trimezieae_corrRes.plot <- ggplot(data = corr_trimezieae, mapping = aes(residuals, total_residuals))+
#  geom_point(position = "jitter")+
#  labs(title = "Trimezieae", subtitle = paste("r =", trimezieae_corrRes))+
#  xlab("Residuals")+
#  ylab("Total residuals")+
#  theme_bw(base_size = 14)+
#  theme(plot.title = element_text(face = "italic"))
#png("plots/synthesis/pd/res_corr/trimezieae.png",
#    height = 4, width = 4, units = 'in', res=300); trimezieae_corrRes.plot; dev.off()

#Marcetieae
marcetieae_corrRes <- as.character(formatC(cor(corr_marcetieae$residuals, 
                                         corr_marcetieae$total_residuals, use = "na.or.complete"))) 
marcetieae_corrRes.plot <- ggplot(data = corr_marcetieae, mapping = aes(residuals, total_residuals))+
  geom_point(position = "jitter")+
  labs(title = "Marcetieae", subtitle = paste("r =", marcetieae_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw(base_size = 14)
png("plots/synthesis/pd/res_corr/marcetieae.png",
    height = 4, width = 4, units = 'in', res=300); marcetieae_corrRes.plot; dev.off()

#=====================================================================================================#
#=====================================================================================================#
#=====================================================================================================#

#==========#
# BOXPLOTS #
#==========#

mimosa_res <- data.frame(group = "Mimosa", res = mimosa_pd$residuals)
diplusodon_res <- data.frame(group = "Diplusodon", res = diplusodon_pd$residuals)
coman_res <- data.frame(group = "Comanthera", res = coman_pd$residuals)
vel_res <- data.frame(group = "Velloziaceae", res = vel_pd$residuals)
lych_res <- data.frame(group = "Lychnophorinae", res = lych_pd$residuals)
paep_res <- data.frame(group = "Paepalanthus", res = paep_pd$residuals)
calli_res <- data.frame(group = "Calliandra", res = calli_pd$residuals)
minaria_res <- data.frame(group = "Minaria", res = minaria_pd$residuals)
marcetieae_res <- data.frame(group = "Marcetieae", res = marcetieae_pd$residuals)
trimezieae_res <- data.frame(group = "Trimezieae", res = trimezieae_pd$residuals)
habenaria_res <- data.frame(group = "Habenaria", res = habenaria_pd$residuals)
catt_res <- data.frame(group = "Cattleya", res = catt_pd$residuals)

residuals <- rbind(mimosa_res, diplusodon_res, coman_res, vel_res, lych_res,
                   paep_res, calli_res, minaria_res, marcetieae_res, trimezieae_res,
                   habenaria_res, catt_res)

library(ggplot2)

ggplot(residuals) +
 aes(x = "", y = res, group = group) +
 geom_boxplot(shape = "circle", fill = "#112446") +
 theme_minimal()

#Violin plot for all groups combined

library(ggplot2)

ggplot(conc) +
 aes(x = "", y = total_residuals) +
 geom_violin(adjust = 0.5, scale = "area", fill = "#E42205") +
 labs(x = "Grupos combinados", y = "Resíduos") +
 theme_bw()


#=====================================================================================================#
#=====================================================================================================#
#=====================================================================================================#

#====================#
# WEIGHTED SYNTHESIS #
#====================#

#Linear model and residuals
set.seed(7)
pd_lm.total <- lm(pd.obs_w ~ ntaxa_w, conc_w)
conc_w$total_residuals <- pd_lm.total$residuals

#Merging synthesis with spatial grids
pdPoly <- merge(grids_cr, conc_w, by.x = "id")

#=======#
# Plots #
#=======#

#SR
sr_plot <- spplot(pdPoly,
                  zcol = "ntaxa_w",
                  xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                  par.settings=list(fontsize = list(text = 18)),
                  colorkey = TRUE, sp.layout = list(list(br, 
                                                         fill = "gray")), 
                  col.regions = rev(magma(16)), scales = list(draw = FALSE))
#png("plots/synthesis/weighted_pd/sr_plot.png",
 #   height = 4, width = 4, units = 'in', res=300); sr_plot; dev.off()

#PD
pd_plot <- spplot(pdPoly,
                  zcol = "pd.obs_w",
                  xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                  par.settings=list(fontsize = list(text = 18)),
                  colorkey = TRUE, sp.layout = list(list(br, 
                                                         fill = "gray")), 
                  col.regions = rev(magma(16)), scales = list(draw = FALSE))
#png("plots/synthesis/weighted_pd/pd_plot.png",
 #   height = 4, width = 4, units = 'in', res=300); pd_plot; dev.off()

#Residuals
pdres_plot <- spplot(pdPoly,
                     zcol = "total_residuals",
                     xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5),
                     at = seq(-0.086, 0.086, length.out = 16),
                     par.settings=list(fontsize = list(text = 18)),
                     colorkey = TRUE, sp.layout = list(list(br, 
                                                            fill = "gray")), 
                     col.regions = viridis(16), scales = list(draw = FALSE))
#png("plots/synthesis/weighted_pd/pdres_plot.png",
 #   height = 4, width = 4, units = 'in', res=300); pdres_plot; dev.off()

#==============#
# CORRELATIONS #
#==============#

#Joining conc_watenated dataset and individual datasets
corr_mimosa <- plyr::join(conc_w[ , c("id", "total_residuals")], mimosa_pd, by = "id")
corr_diplusodon <- plyr::join(conc_w[ , c("id", "total_residuals")], diplusodon_pd, by = "id")
corr_coman <- plyr::join(conc_w[ , c("id", "total_residuals")], coman_pd, by = "id")
corr_vel <- plyr::join(conc_w[ , c("id", "total_residuals")], vel_pd, by = "id")
corr_lych <- plyr::join(conc_w[ , c("id", "total_residuals")], lych_pd, by = "id")
corr_paep <- plyr::join(conc_w[ , c("id", "total_residuals")], paep_pd, by = "id")

#=======#
# Plots #
#=======#

#Mimosa
mimosa_corrRes <- as.character(formatC(cor(corr_mimosa$residuals, 
                                           corr_mimosa$total_residuals, use = "na.or.complete"))) 
mimosa_corrRes.plot <- ggplot(data = corr_mimosa, mapping = aes(residuals, total_residuals))+
  geom_point()+
  labs(title = "Mimosa", subtitle = paste("r =", mimosa_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw()+
  theme(plot.title = element_text(face = "italic"))
png("plots/synthesis/weighted_pd/res_corr/mimosa.png",
    height = 4, width = 4, units = 'in', res=300); mimosa_corrRes.plot; dev.off()

#Diplusodon
diplusodon_corrRes <- as.character(formatC(cor(corr_diplusodon$residuals, 
                                               corr_diplusodon$total_residuals, use = "na.or.complete"))) 
diplusodon_corrRes.plot <- ggplot(data = corr_diplusodon, mapping = aes(residuals, total_residuals))+
  geom_point()+
  labs(title = "Diplusodon", subtitle = paste("r =", diplusodon_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw()+
  theme(plot.title = element_text(face = "italic"))
png("plots/synthesis/weighted_pd/res_corr/diplusodon.png",
    height = 4, width = 4, units = 'in', res=300); diplusodon_corrRes.plot; dev.off()

#Comanthera
coman_corrRes <- as.character(formatC(cor(corr_coman$residuals, 
                                          corr_coman$total_residuals, use = "na.or.complete"))) 
coman_corrRes.plot <- ggplot(data = corr_coman, mapping = aes(residuals, total_residuals))+
  geom_point()+
  labs(title = "Comanthera", subtitle = paste("r =", coman_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw()+
  theme(plot.title = element_text(face = "italic"))
png("plots/synthesis/weighted_pd/res_corr/coman.png",
    height = 4, width = 4, units = 'in', res=300); coman_corrRes.plot; dev.off()

#Velloziaceae
vel_corrRes <- as.character(formatC(cor(corr_vel$residuals, 
                                        corr_vel$total_residuals, use = "na.or.complete"))) 
vel_corrRes.plot <- ggplot(data = corr_vel, mapping = aes(residuals, total_residuals))+
  geom_point()+
  labs(title = "Velloziaceae", subtitle = paste("r =", vel_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw()+
  theme()
png("plots/synthesis/weighted_pd/res_corr/vel.png",
    height = 4, width = 4, units = 'in', res=300); vel_corrRes.plot; dev.off()

#Lychnophorinae
lych_corrRes <- as.character(formatC(cor(corr_lych$residuals, 
                                         corr_lych$total_residuals, use = "na.or.complete"))) 
lych_corrRes.plot <- ggplot(data = corr_lych, mapping = aes(residuals, total_residuals))+
  geom_point()+
  labs(title = "Lychnophorinae", subtitle = paste("r =", lych_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw()+
  theme()
png("plots/synthesis/weighted_pd/res_corr/lych.png",
    height = 4, width = 4, units = 'in', res=300); lych_corrRes.plot; dev.off()

#Paepalanthus
paep_corrRes <- as.character(formatC(cor(corr_paep$residuals, 
                                         corr_paep$total_residuals, use = "na.or.complete"))) 
paep_corrRes.plot <- ggplot(data = corr_paep, mapping = aes(residuals, total_residuals))+
  geom_point()+
  labs(title = "Paepalanthus", subtitle = paste("r =", paep_corrRes))+
  xlab("Residuals")+
  ylab("Total residuals")+
  theme_bw()+
  theme(plot.title = element_text(face = "italic"))
png("plots/synthesis/weighted_pd/res_corr/paep.png",
    height = 4, width = 4, units = 'in', res=300); paep_corrRes.plot; dev.off()