library(iNEXT)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

calli_cr <- read.csv("datasets/Calliandra/calli_cr.csv")

# abundance matrix
calli_matrix <- matrix(data = NA, nrow = length(unique(calli_cr$id_grid)), 
                        ncol = length(unique(calli_cr$gen_sp)))
calli_matrix <- as.data.frame(calli_matrix)
colnames(calli_matrix) <- unique(calli_cr$gen_sp)
rownames(calli_matrix) <- unique(calli_cr$id_grid)
for(i in 1:nrow(calli_matrix)){
  for(j in 1:ncol(calli_matrix)){
    if(colnames(calli_matrix)[j] %in% calli_cr$gen_sp[calli_cr$id_grid == rownames(calli_matrix)[i]]){
      calli_matrix[i, j] <- nrow(calli_cr[calli_cr$gen_sp == colnames(calli_matrix)[j] & calli_cr$id_grid == as.numeric(rownames(calli_matrix)[i]), ])
    } else {
      calli_matrix[i, j] <- 0  
    }
  }
}

# Transposing matrix
calli_matrix.t <- t(calli_matrix)

# Removing cells with only one species
calli_matrix.t <- calli_matrix.t[ , names(which(colSums(calli_matrix.t == 0) < nrow(calli_matrix.t) - 1))]

# Defining sample sizes
m <- c(1, 2, 10, 20, 50, 100, 200, 500, 1000, 2000)

# Running iNEXT
i.next <- iNEXT(x = calli_matrix.t, datatype = "abundance", size = m)

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
  labs(title = "Calliandra", subtitle = paste("r =", corrSr.ObsEst, "\np < 0.05, R² = 0.94, slope = 0.75"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 21)+
  theme(plot.title = element_text(face = "italic"))+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Calliandra/calli-corrSrObsEst_plot.pdf"); corrSr.ObsEst_plot; dev.off()
