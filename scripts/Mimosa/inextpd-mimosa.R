library(iNextPD)
library(ape)
library(ade4)
library(TreeTools)
library(reshape2)
library(tidyverse)

# Set wd
setwd("B:/yagob/GoogleDrive/Academia/R-MSc")

# Reading data set
mimosa_cr <- read.csv("datasets/Mimosa/mimosa_cr.csv")

# Reading tree
tree <- read.nexus("trees/Mimosa/pruned_tree-mimosa.nex")
mimosa.phy <- NewickTree(tree)
mimosa.phy <- gsub(" ", "_", mimosa.phy)
mimosa.phy <- ade4::newick2phylog(mimosa.phy)

# Removing species not sampled in the tree
mimosa_cr <- mimosa_cr[which(mimosa_cr$gen_sp %in% tree$tip.label), ]

# abundance matrix
mimosa_matrix <- matrix(data = NA, nrow = length(unique(mimosa_cr$id_grid)), 
                        ncol = length(unique(mimosa_cr$gen_sp)))
mimosa_matrix <- as.data.frame(mimosa_matrix)
colnames(mimosa_matrix) <- unique(mimosa_cr$gen_sp)
rownames(mimosa_matrix) <- unique(mimosa_cr$id_grid)
for(i in 1:nrow(mimosa_matrix)){
  for(j in 1:ncol(mimosa_matrix)){
    if(colnames(mimosa_matrix)[j] %in% mimosa_cr$gen_sp[mimosa_cr$id_grid == rownames(mimosa_matrix)[i]]){
      mimosa_matrix[i, j] <- nrow(mimosa_cr[mimosa_cr$gen_sp == colnames(mimosa_matrix)[j] & mimosa_cr$id_grid == as.numeric(rownames(mimosa_matrix)[i]), ])
    } else {
      mimosa_matrix[i, j] <- 0  
    }
  }
}

# Transposing matrix
mimosa_matrix.t <- t(mimosa_matrix)

# Removing cells with only one species
mimosa_matrix.t <- mimosa_matrix.t[ , names(which(colSums(mimosa_matrix.t == 0) < nrow(mimosa_matrix.t) - 1))]

# Running iNextPD and preparing data frame
i.nextpd <- iNextPD(x = mimosa_matrix.t, labels = row.names(mimosa_matrix.t), phy = mimosa.phy)
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
  labs(title = "Mimosa", subtitle = paste("r =", corrPd.ObsEst, "\np < 0.05, R² = 0.97"))+
  xlab("Estimated")+
  ylab("Observed")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(face = "italic"))+
  geom_smooth(method='lm', formula= y~x)
cairo_pdf("figures/Mimosa/mimosa-corrPdObsEst_plot.pdf"); corrPd.ObsEst_plot; dev.off()

