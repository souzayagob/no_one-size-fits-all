library(iNEXT)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

diplusodon_cr <- read.csv("datasets/Diplusodon/diplusodon_cr.csv")

# abundance matrix
diplusodon_matrix <- matrix(data = NA, nrow = length(unique(diplusodon_cr$id_grid)), 
                        ncol = length(unique(diplusodon_cr$gen_sp)))
diplusodon_matrix <- as.data.frame(diplusodon_matrix)
colnames(diplusodon_matrix) <- unique(diplusodon_cr$gen_sp)
rownames(diplusodon_matrix) <- unique(diplusodon_cr$id_grid)
for(i in 1:nrow(diplusodon_matrix)){
  for(j in 1:ncol(diplusodon_matrix)){
    if(colnames(diplusodon_matrix)[j] %in% diplusodon_cr$gen_sp[diplusodon_cr$id_grid == rownames(diplusodon_matrix)[i]]){
      diplusodon_matrix[i, j] <- nrow(diplusodon_cr[diplusodon_cr$gen_sp == colnames(diplusodon_matrix)[j] & diplusodon_cr$id_grid == as.numeric(rownames(diplusodon_matrix)[i]), ])
    } else {
      diplusodon_matrix[i, j] <- 0  
    }
  }
}

# Transposing matrix
diplusodon_matrix.t <- t(diplusodon_matrix)

# Removing cells with only one species
diplusodon_matrix.t <- diplusodon_matrix.t[ , names(which(colSums(diplusodon_matrix.t == 0) < nrow(diplusodon_matrix.t) - 1))]

# Defining sample sizes
m <- c(1, 2, 10, 20, 50, 100, 200, 500, 1000, 2000)

# Running iNEXT
i.next <- iNEXT(x = diplusodon_matrix.t, datatype = "abundance", size = m)

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
  labs(title = "Diplusodon", subtitle = paste("r =", corrSr.ObsEst, "\np < 0.05, R² = 0.93, slope = 0.74"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 21)+
  theme(plot.title = element_text(face = "italic"))+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Diplusodon/diplusodon-corrSrObsEst_plot.pdf"); corrSr.ObsEst_plot; dev.off()
