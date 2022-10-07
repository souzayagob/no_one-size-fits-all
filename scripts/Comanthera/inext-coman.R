library(iNEXT)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

coman_cr <- read.csv("datasets/Comanthera/coman_cr.csv")

# abundance matrix
coman_matrix <- matrix(data = NA, nrow = length(unique(coman_cr$id_grid)), 
                        ncol = length(unique(coman_cr$gen_sp)))
coman_matrix <- as.data.frame(coman_matrix)
colnames(coman_matrix) <- unique(coman_cr$gen_sp)
rownames(coman_matrix) <- unique(coman_cr$id_grid)
for(i in 1:nrow(coman_matrix)){
  for(j in 1:ncol(coman_matrix)){
    if(colnames(coman_matrix)[j] %in% coman_cr$gen_sp[coman_cr$id_grid == rownames(coman_matrix)[i]]){
      coman_matrix[i, j] <- nrow(coman_cr[coman_cr$gen_sp == colnames(coman_matrix)[j] & coman_cr$id_grid == as.numeric(rownames(coman_matrix)[i]), ])
    } else {
      coman_matrix[i, j] <- 0  
    }
  }
}

# Transposing matrix
coman_matrix.t <- t(coman_matrix)

# Removing cells with only one species
coman_matrix.t <- coman_matrix.t[ , names(which(colSums(coman_matrix.t == 0) < nrow(coman_matrix.t) - 1))]

# Defining sample sizes
m <- c(1, 2, 10, 20, 50, 100, 200, 500, 1000, 2000)

# Running iNEXT
i.next <- iNEXT(x = coman_matrix.t, datatype = "abundance", size = m)

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
  labs(title = "Comanthera", subtitle = paste("r =", corrSr.ObsEst, "\np < 0.05, R² = 0.85, slope = 0.70"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 21)+
  theme(plot.title = element_text(face = "italic"))+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Comanthera/coman-corrSrObsEst_plot.pdf"); corrSr.ObsEst_plot; dev.off()
