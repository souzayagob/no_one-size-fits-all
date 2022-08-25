library(iNEXT)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

paep_cr <- read.csv("datasets/Paepalanthus/paep_cr.csv")

# abundance matrix
paep_matrix <- matrix(data = NA, nrow = length(unique(paep_cr$id_grid)), 
                        ncol = length(unique(paep_cr$gen_sp)))
paep_matrix <- as.data.frame(paep_matrix)
colnames(paep_matrix) <- unique(paep_cr$gen_sp)
rownames(paep_matrix) <- unique(paep_cr$id_grid)
for(i in 1:nrow(paep_matrix)){
  for(j in 1:ncol(paep_matrix)){
    if(colnames(paep_matrix)[j] %in% paep_cr$gen_sp[paep_cr$id_grid == rownames(paep_matrix)[i]]){
      paep_matrix[i, j] <- nrow(paep_cr[paep_cr$gen_sp == colnames(paep_matrix)[j] & paep_cr$id_grid == as.numeric(rownames(paep_matrix)[i]), ])
    } else {
      paep_matrix[i, j] <- 0  
    }
  }
}

# Transposing matrix
paep_matrix.t <- t(paep_matrix)

# Removing cells with only one species
paep_matrix.t <- paep_matrix.t[ , names(which(colSums(paep_matrix.t == 0) < nrow(paep_matrix.t) - 1))]

# Defining sample sizes
m <- c(1, 2, 10, 20, 50, 100, 200, 500, 1000, 2000)

# Running iNEXT
i.next <- iNEXT(x = paep_matrix.t, datatype = "abundance", size = m)

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
  labs(title = "Paepalanthus", subtitle = paste("r =", corrSr.ObsEst, "\np < 0.05, R² = 0.99"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Paepalanthus/paep-corrSrObsEst_plot.pdf"); corrSr.ObsEst_plot; dev.off()