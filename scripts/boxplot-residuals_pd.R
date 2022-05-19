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

#=============#
# Data frames #
#=============#

#Data frame with all grid ids
grids_df <- as.data.frame(grids_cr@data)

#Data frame containing sum of SR and PD for each ID. Weighted (conc_w) nd non-weigthted (conc)
conc <- plyr::join(grids_df, mimosa_pd, by = "id")
conc <- plyr::join(conc[ , -which(colnames(conc) == "layer")], diplusodon_pd, by = "id")
conc <- plyr::join(conc, coman_pd, by = "id")
conc <- plyr::join(conc, vel_pd, by = "id")
conc <- plyr::join(conc, lych_pd, by = "id")
conc <- plyr::join(conc, paep_pd, by = "id")
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

#=====================================================================================================#

#=========#
# BOXPLOT #
#=========#

#Linear model and residuals - Mimosa
set.seed(7)
pd_lm.mimosa <- lm(pd.obs ~ ntaxa, mimosa_pd)
mimosa_pd$residuals <- pd_lm.mimosa$residuals

#Linear model and residuals - Diplusodon
set.seed(7)
pd_lm.diplusodon <- lm(pd.obs ~ ntaxa, diplusodon_pd)
diplusodon_pd$residuals <- pd_lm.diplusodon$residuals

#Linear model and residuals - Comanthera
set.seed(7)
pd_lm.coman <- lm(pd.obs ~ ntaxa, coman_pd)
coman_pd$residuals <- pd_lm.coman$residuals

#Linear model and residuals _ Velloziaceae
set.seed(7)
pd_lm.vel <- lm(pd.obs ~ ntaxa, vel_pd)
vel_pd$residuals <- pd_lm.vel$residuals

#Linear model and residuals - Lychnophorinae
set.seed(7)
pd_lm.lych <- lm(pd.obs ~ ntaxa, lych_pd)
lych_pd$residuals <- pd_lm.lych$residuals

#Linear model and residuals - Paepalanthus
set.seed(7)
pd_lm.paep <- lm(pd.obs ~ ntaxa, paep_pd)
paep_pd$residuals <- pd_lm.paep$residuals

#Boxplots
box_df <- data.frame("residuals" = mimosa_pd$residuals, "group" = "Mimosa")
box_df <- rbind(box_df, data.frame("residuals" = diplusodon_pd$residuals, "group" = "Diplusodon"))
box_df <- rbind(box_df, data.frame("residuals" = coman_pd$residuals, "group" = "Comanthera"))
box_df <- rbind(box_df, data.frame("residuals" = vel_pd$residuals, "group" = "Velloziaceae"))
box_df <- rbind(box_df, data.frame("residuals" = lych_pd$residuals, "group" = "Lychnophorinae"))
box_df <- rbind(box_df, data.frame("residuals" = paep_pd$residuals, "group" = "Paepalanthus"))
box_df <- rbind(box_df, data.frame("residuals" = conc$total_residuals, "group" = "All groups"))

boxplot(residuals ~ group, box_df)