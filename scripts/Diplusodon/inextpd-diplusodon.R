library(iNextPD)
library(ape)
library(ade4)
library(TreeTools)
library(reshape2)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

# Reading data set
diplusodon_cr <- read.csv("datasets/Diplusodon/diplusodon_cr.csv")

# Reading tree
tree <- read.nexus("trees/Diplusodon/pruned_tree-diplusodon.nex")
diplusodon.phy <- NewickTree(tree) 
diplusodon.phy <- gsub(" ", "_", diplusodon.phy)
diplusodon.phy <- ade4::newick2phylog(diplusodon.phy)

# Removing species not sampled in the tree
diplusodon_cr <- diplusodon_cr[which(diplusodon_cr$gen_sp %in% tree$tip.label), ]

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

# Running iNextPD and preparing data frame
i.nextpd <- iNextPD(x = diplusodon_matrix.t, labels = row.names(diplusodon_matrix.t), phy = diplusodon.phy)
masy.est <- i.nextpd$AsyPDEst
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
  labs(title = "Diplusodon", subtitle = paste("r =", corrPd.ObsEst, "\np < 0.05, R² = 0.97"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Diplusodon/diplusodon-corrPdObsEst_plot.pdf"); corrPd.ObsEst_plot; dev.off()