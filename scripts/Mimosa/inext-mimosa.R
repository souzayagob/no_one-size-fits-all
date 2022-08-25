library(iNEXT)

mimosa_cr <- read.csv("datasets/Mimosa/mimosa_cr.csv")

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

# Defining sample sizes
m <- c(1, 2, 10, 20, 50, 100, 200, 500, 1000, 2000)

# Running iNEXT
i.next <- iNEXT(x = mimosa_matrix.t, datatype = "abundance", size = m)

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