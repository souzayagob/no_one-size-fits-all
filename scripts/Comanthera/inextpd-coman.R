library(iNextPD)
library(ape)
library(ade4)
library(TreeTools)
library(reshape2)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

# Reading data set
coman_cr <- read.csv("datasets/Comanthera/coman_cr.csv")

# Reading tree
tree <- read.nexus("trees/Comanthera/pruned_tree-coman.nex")
coman.phy <- NewickTree(tree)
coman.phy <- gsub(" ", "_", coman.phy)
coman.phy <- ade4::newick2phylog(coman.phy)

# Removing species not sampled in the tree
coman_cr <- coman_cr[which(coman_cr$gen_sp %in% tree$tip.label), ]

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

# Running iNextPD and preparing data frame
i.nextpd <- iNextPD(x = as.data.frame(coman_matrix.t), labels = row.names(coman_matrix.t), phy = coman.phy,
                    datatype = "abundance")
asy.est <- i.nextpd$AsyPDEst
asy.est <- as.data.frame(asy.est)
asy.est <- dcast(asy.est, Var1+Var2~Var3, value.var = "Freq")
asy.est <- asy.est[asy.est$Var2 == "q = 0", ]

# Are the observed values correlated with the estimated values? 
pd.lm <- lm(asy.est$Observed ~ asy.est$Estimator)
summary(pd.lm)

# Computing Pearson correlation
corrPd.ObsEst <- as.character(formatC(cor(asy.est$Observed, 
                                          asy.est$Estimator, use = "na.or.complete")))

# Plotting
corrPd.ObsEst_plot <- ggplot(data = asy.est, mapping = aes(Estimator, Observed))+
  geom_jitter()+
  labs(title = "Comanthera", subtitle = paste("r =", corrPd.ObsEst, "\np < 0.05, R² = 0.97"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Comanthera/coman-corrPdObsEst_plot.pdf"); corrPd.ObsEst_plot; dev.off()
