library(iNEXT)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

lych_cr <- read.csv("datasets/Lychnophorinae/lych_cr.csv")

# abundance matrix
lych_matrix <- matrix(data = NA, nrow = length(unique(lych_cr$id_grid)), 
                        ncol = length(unique(lych_cr$gen_sp)))
lych_matrix <- as.data.frame(lych_matrix)
colnames(lych_matrix) <- unique(lych_cr$gen_sp)
rownames(lych_matrix) <- unique(lych_cr$id_grid)
for(i in 1:nrow(lych_matrix)){
  for(j in 1:ncol(lych_matrix)){
    if(colnames(lych_matrix)[j] %in% lych_cr$gen_sp[lych_cr$id_grid == rownames(lych_matrix)[i]]){
      lych_matrix[i, j] <- nrow(lych_cr[lych_cr$gen_sp == colnames(lych_matrix)[j] & lych_cr$id_grid == as.numeric(rownames(lych_matrix)[i]), ])
    } else {
      lych_matrix[i, j] <- 0  
    }
  }
}

# Transposing matrix
lych_matrix.t <- t(lych_matrix)

# Removing cells with only one species
lych_matrix.t <- lych_matrix.t[ , names(which(colSums(lych_matrix.t == 0) < nrow(lych_matrix.t) - 1))]

# Defining sample sizes
m <- c(1, 2, 10, 20, 50, 100, 200, 500, 1000, 2000)

# Running iNEXT
i.next <- iNEXT(x = lych_matrix.t, datatype = "abundance", size = m)

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
  labs(title = "Lychnophorinae", subtitle = paste("r =", corrSr.ObsEst, "\np < 0.05, R² = 0.95"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text())+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Lychnophorinae/lych-corrSrObsEst_plot.pdf"); corrSr.ObsEst_plot; dev.off()
