library(iNEXT)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

vel_cr <- read.csv("datasets/Velloziaceae/vel_cr.csv")

# abundance matrix
vel_matrix <- matrix(data = NA, nrow = length(unique(vel_cr$id_grid)), 
                        ncol = length(unique(vel_cr$gen_sp)))
vel_matrix <- as.data.frame(vel_matrix)
colnames(vel_matrix) <- unique(vel_cr$gen_sp)
rownames(vel_matrix) <- unique(vel_cr$id_grid)
for(i in 1:nrow(vel_matrix)){
  for(j in 1:ncol(vel_matrix)){
    if(colnames(vel_matrix)[j] %in% vel_cr$gen_sp[vel_cr$id_grid == rownames(vel_matrix)[i]]){
      vel_matrix[i, j] <- nrow(vel_cr[vel_cr$gen_sp == colnames(vel_matrix)[j] & vel_cr$id_grid == as.numeric(rownames(vel_matrix)[i]), ])
    } else {
      vel_matrix[i, j] <- 0  
    }
  }
}

# Transposing matrix
vel_matrix.t <- t(vel_matrix)

# Removing cells with only one species
vel_matrix.t <- vel_matrix.t[ , names(which(colSums(vel_matrix.t == 0) < nrow(vel_matrix.t) - 1))]

# Defining sample sizes
m <- c(1, 2, 10, 20, 50, 100, 200, 500, 1000, 2000)

# Running iNEXT
i.next <- iNEXT(x = vel_matrix.t, datatype = "abundance", size = m)

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
  labs(title = "Velloziaceae", subtitle = paste("r =", corrSr.ObsEst, "\np < 0.05, R² = 0.96"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text())+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Velloziaceae/vel-corrSrObsEst_plot.pdf"); corrSr.ObsEst_plot; dev.off()
