library(iNEXT)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

marcetieae_cr <- read.csv("datasets/Marcetieae/marcetieae_cr.csv")

# abundance matrix
marcetieae_matrix <- matrix(data = NA, nrow = length(unique(marcetieae_cr$id_grid)), 
                        ncol = length(unique(marcetieae_cr$gen_sp)))
marcetieae_matrix <- as.data.frame(marcetieae_matrix)
colnames(marcetieae_matrix) <- unique(marcetieae_cr$gen_sp)
rownames(marcetieae_matrix) <- unique(marcetieae_cr$id_grid)
for(i in 1:nrow(marcetieae_matrix)){
  for(j in 1:ncol(marcetieae_matrix)){
    if(colnames(marcetieae_matrix)[j] %in% marcetieae_cr$gen_sp[marcetieae_cr$id_grid == rownames(marcetieae_matrix)[i]]){
      marcetieae_matrix[i, j] <- nrow(marcetieae_cr[marcetieae_cr$gen_sp == colnames(marcetieae_matrix)[j] & marcetieae_cr$id_grid == as.numeric(rownames(marcetieae_matrix)[i]), ])
    } else {
      marcetieae_matrix[i, j] <- 0  
    }
  }
}

# Transposing matrix
marcetieae_matrix.t <- t(marcetieae_matrix)

# Removing cells with only one species
marcetieae_matrix.t <- marcetieae_matrix.t[ , names(which(colSums(marcetieae_matrix.t == 0) < nrow(marcetieae_matrix.t) - 1))]

# Defining sample sizes
m <- c(1, 2, 10, 20, 50, 100, 200, 500, 1000, 2000)

# Running iNEXT
i.next <- iNEXT(x = marcetieae_matrix.t, datatype = "abundance", size = m)

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
  labs(title = "Marcetieae", subtitle = paste("r =", corrSr.ObsEst, "\np < 0.05, R² = 0.82"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text())+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Marcetieae/marcetieae-corrSrObsEst_plot.pdf"); corrSr.ObsEst_plot; dev.off()

# How the correlation and linear model change when we remove the outlier?
no.out <- asy.sr[asy.sr$Site != "A", ]
no.out_lm <- lm(no.out$Observed ~ no.out$Estimator)
summary(no.out_lm)
corrNoOut <- as.character(formatC(cor(no.out$Observed, 
                                          no.out$Estimator, use = "na.or.complete")))
corrNoOut_plot <- ggplot(data = no.out, mapping = aes(Estimator, Observed))+
  geom_jitter()+
  labs(title = "Marcetieae", subtitle = paste("r =", corrNoOut, "\np < 0.05, R² = 0.97"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text())+
  geom_smooth(method='lm', formula= y~x)
