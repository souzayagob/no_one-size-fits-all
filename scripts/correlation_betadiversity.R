#=====================================================================================================#

#=======#
# INPUT #
#=======#

#========#
# Mimosa #
#========#

#Jaccard distance
mimosa_j <- as.matrix(read.csv("results/Mimosa/jaccard_matrix.csv", row.names = 1))
colnames(mimosa_j) <- gsub("X", "", colnames(mimosa_j))

#UniFrac
mimosa_u <- as.matrix(read.csv("results/Mimosa/unifrac_matrix.csv", row.names = 1))
colnames(mimosa_u) <- gsub("X", "", colnames(mimosa_u))

#============#
# Diplusodon #
#============#

#Jaccard distance
diplusodon_j <- as.matrix(read.csv("results/Diplusodon/jaccard_matrix.csv", row.names = 1))
colnames(diplusodon_j) <- gsub("X", "", colnames(diplusodon_j))

#UniFrac
diplusodon_u <- as.matrix(read.csv("results/Diplusodon/unifrac_matrix.csv", row.names = 1))
colnames(diplusodon_u) <- gsub("X", "", colnames(diplusodon_u))

#============#
# Comanthera #
#============#

#Jaccard distance
coman_j <- as.matrix(read.csv("results/Comanthera/jaccard_matrix.csv", row.names = 1))
colnames(coman_j) <- gsub("X", "", colnames(coman_j))

#UniFrac
coman_u <- as.matrix(read.csv("results/Comanthera/unifrac_matrix.csv", row.names = 1))
colnames(coman_u) <- gsub("X", "", colnames(coman_u))

#==============#
# Velloziaceae #
#==============#

#Jaccard distance
vel_j <- as.matrix(read.csv("results/Velloziaceae/jaccard_matrix.csv", row.names = 1))
colnames(vel_j) <- gsub("X", "", colnames(vel_j))

#UniFrac
vel_u <- as.matrix(read.csv("results/Velloziaceae/unifrac_matrix.csv", row.names = 1))
colnames(vel_u) <- gsub("X", "", colnames(vel_u))

#================#
# Lychnophorinae #
#================#

#Jaccard distance
lych_j <- as.matrix(read.csv("results/Lychnophorinae/jaccard_matrix.csv", row.names = 1))
colnames(lych_j) <- gsub("X", "", colnames(lych_j))

#UniFrac
lych_u <- as.matrix(read.csv("results/Lychnophorinae/unifrac_matrix.csv", row.names = 1))
colnames(lych_u) <- gsub("X", "", colnames(lych_u))

#==============#
# Paepalanthus #
#==============#

#Jaccard distance
paep_j <- as.matrix(read.csv("results/Paepalanthus/jaccard_matrix.csv", row.names = 1))
colnames(paep_j) <- gsub("X", "", colnames(paep_j))

#UniFrac
paep_u <- as.matrix(read.csv("results/Paepalanthus/unifrac_matrix.csv", row.names = 1))
colnames(paep_u) <- gsub("X", "", colnames(paep_u))

#============#
# All groups #
#============#

all_j <- as.matrix(read.csv("results/synthesis/jaccard_matrix.csv", row.names = 1))
colnames(all_j) <- gsub("X", "", colnames(all_j))


#=====================================================================================================#



#==============#
# Correlations #
#==============#

#======================#
# Jaccard ~ all.Jacard #
#======================#

#Mimosa
cor(c(mimosa_j), c(all_j[which(rownames(all_j) %in% rownames(mimosa_j))
                         , which(colnames(all_j) %in% colnames(mimosa_j))]))

#Diplusodon
cor(c(diplusodon_j), c(all_j[which(rownames(all_j) %in% rownames(diplusodon_j))
                         , which(colnames(all_j) %in% colnames(diplusodon_j))]))

#Comanthera
cor(c(coman_j), c(all_j[which(rownames(all_j) %in% rownames(coman_j))
                         , which(colnames(all_j) %in% colnames(coman_j))]))

#Velloziaceae
cor(c(vel_j), c(all_j[which(rownames(all_j) %in% rownames(vel_j))
                         , which(colnames(all_j) %in% colnames(vel_j))]))

#Lychnophorinae
cor(c(lych_j), c(all_j[which(rownames(all_j) %in% rownames(lych_j))
                         , which(colnames(all_j) %in% colnames(lych_j))]))

#Paepalanthus
cor(c(paep_j), c(all_j[which(rownames(all_j) %in% rownames(paep_j))
                         , which(colnames(all_j) %in% colnames(paep_j))]))

#===================#
# Jaccard - Unifrac #
#===================#

#Mimosa
cor(c(mimosa_j), c(mimosa_u))

#Diplusodon
cor(c(diplusodon_j), c(diplusodon_u))

#Comanthera
cor(c(coman_j), c(coman_u))

#Velloziaceae
cor(c(vel_j), c(vel_u))

#Lychnophorinae
cor(c(lych_j), c(lych_u))

#Paepalanthus
cor(c(paep_j), c(paep_u))



