library(iNEXT)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

# Mimosa
mimosa <- read.csv("datasets/Mimosa/mimosa_cr.csv")

# Diplusodon
diplusodon <- read.csv("datasets/Diplusodon/diplusodon_cr.csv")

# Comanthera
coman <- read.csv("datasets/Comanthera/coman_cr.csv")

# Velloziaceae
vel <- read.csv("datasets/Velloziaceae/vel_cr.csv")

# Lychnophorinae
lych <- read.csv("datasets/Lychnophorinae/lych_cr.csv")

# Paepalanthus
paep <- read.csv("datasets/Paepalanthus/paep_cr.csv")

# Calliandra
calli <- read.csv("datasets/Calliandra/calli_cr.csv")

# Minaria
minaria <- read.csv("datasets/Minaria/minaria_cr.csv")

# Marcetieae
marcetieae <- read.csv("datasets/Marcetieae/marcetieae_cr.csv")

# Concatenating
all_groups <- rbind(mimosa,
                    diplusodon,
                    coman,
                    vel,
                    lych,
                    paep,
                    calli,
                    minaria,
                    marcetieae)

# Abundance matrix
all_matrix <- matrix(data = NA, nrow = length(unique(all_groups$id_grid)), 
                     ncol = length(unique(all_groups$gen_sp)))
all_matrix <- as.data.frame(all_matrix)
colnames(all_matrix) <- unique(all_groups$gen_sp)
rownames(all_matrix) <- unique(all_groups$id_grid)
for(i in 1:nrow(all_matrix)){
  for(j in 1:ncol(all_matrix)){
    if(colnames(all_matrix)[j] %in% all_groups$gen_sp[all_groups$id_grid == rownames(all_matrix)[i]]){
      all_matrix[i, j] <- nrow(all_groups[all_groups$gen_sp == colnames(all_matrix)[j] & all_groups$id_grid == as.numeric(rownames(all_matrix)[i]), ])
    } else {
      all_matrix[i, j] <- 0  
    }
  }
}

# Transposing matrix
all_matrix.t <- t(all_matrix)

# Removing cells with only one species
all_matrix.t <- all_matrix.t[ , names(which(colSums(all_matrix.t == 0) < nrow(all_matrix.t) - 1))]

# Defining sample sizes
m <- c(1, 2, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000)

# Running iNEXT
i.next <- iNEXT(x = all_matrix.t, datatype = "abundance", size = m)

# Checking iNextEst
i.next$iNextEst

# Checking DataInfo
i.next$DataInfo

# Data frame with asymptotic estimations for diversity metrics
asy <- i.next$AsyEst

# Data frame with asymptotic estimations for species richness
asy.sr <- asy[asy$Diversity == "Species richness", ]

# Are the observed values correlated with the estimated values? 
sr.lm <- lm(asy.sr$Observed ~ asy.sr$Estimator)
summary(sr.lm)

# Computing Pearson correlation
corrSr.ObsEst <- as.character(formatC(cor(asy.sr$Observed, 
                                          asy.sr$Estimator, use = "na.or.complete")))

# Plotting
corrSr.ObsEst_plot <- ggplot(data = asy.sr, mapping = aes(Estimator, Observed))+
  geom_jitter()+
  labs(title = "All groups", subtitle = paste("r =", corrSr.ObsEst, "\np < 0.05, R² = 0.99, slope = 0.92"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text())+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/all_groups-corrSrObsEst_plot.pdf"); corrSr.ObsEst_plot; dev.off()
