library(iNEXT)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

minaria_cr <- read.csv("datasets/Minaria/minaria_cr.csv")

# abundance matrix
minaria_matrix <- matrix(data = NA, nrow = length(unique(minaria_cr$id_grid)), 
                        ncol = length(unique(minaria_cr$gen_sp)))
minaria_matrix <- as.data.frame(minaria_matrix)
colnames(minaria_matrix) <- unique(minaria_cr$gen_sp)
rownames(minaria_matrix) <- unique(minaria_cr$id_grid)
for(i in 1:nrow(minaria_matrix)){
  for(j in 1:ncol(minaria_matrix)){
    if(colnames(minaria_matrix)[j] %in% minaria_cr$gen_sp[minaria_cr$id_grid == rownames(minaria_matrix)[i]]){
      minaria_matrix[i, j] <- nrow(minaria_cr[minaria_cr$gen_sp == colnames(minaria_matrix)[j] & minaria_cr$id_grid == as.numeric(rownames(minaria_matrix)[i]), ])
    } else {
      minaria_matrix[i, j] <- 0  
    }
  }
}

# Transposing matrix
minaria_matrix.t <- t(minaria_matrix)

# Removing cells with only one species
minaria_matrix.t <- minaria_matrix.t[ , names(which(colSums(minaria_matrix.t == 0) < nrow(minaria_matrix.t) - 1))]

# Defining sample sizes
m <- c(1, 2, 10, 20, 50, 100, 200, 500, 1000, 2000)

# Running iNEXT
i.next <- iNEXT(x = minaria_matrix.t, datatype = "abundance", size = m)

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
  labs(title = "Minaria", subtitle = paste("r =", corrSr.ObsEst, "\np < 0.05, R² = 0.98, slope = 0.86"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 21)+
  theme(plot.title = element_text(face = "italic"))+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Minaria/minaria-corrSrObsEst_plot.pdf"); corrSr.ObsEst_plot; dev.off()
