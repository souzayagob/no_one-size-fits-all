library(tidyverse) 
library(data.table) 
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 

source('C:/Users/Yago/Dropbox/Mestrado/Bancos_organizing/Scripts/functions.R') #Windows
source('~/Dropbox/Mestrado/R-MSc/scripts/functions.R') #Ubuntu
#======================================================================================#

#==================#
# READING DATASETS #
#==================#

#======#
# GBIF #
#======#

#Reading gbif data
paep_gbif <- fread(file = "datasets/Paepalanthus/0046916-200613084148143-gbif_paepalanthus/occurrence.txt",
                   na.strings = c("", NA), stringsAsFactors = FALSE, encoding = "UTF-8")

#Selecting important attributes
paep_gbif <- paep_gbif %>% dplyr::select(institutionCode,
                                  collectionCode,
                                  catalogNumber,
                                  genus,
                                  specificEpithet,
                                  infraspecificEpithet,
                                  basisOfRecord,
                                  identifiedBy,
                                  recordNumber,
                                  recordedBy,
                                  stateProvince,
                                  county,
                                  municipality,
                                  locality,
                                  decimalLongitude,
                                  decimalLatitude)

#Renaming attributes
paep_gbif <- paep_gbif %>% rename("species" = specificEpithet,
                                  "institutioncode" = institutionCode,
                                  "collectioncode" = collectionCode,
                                  "catalognumber" = catalogNumber,
                                  "basisofrecord" = basisOfRecord,
                                  "identifiedby" = identifiedBy,
                                  "collector" = recordedBy,
                                  "collectornumber" = recordNumber,
                                  "stateprovince" = stateProvince,
                                  "longitude" = decimalLongitude,
                                  "latitude" = decimalLatitude,
                                  "subspecies" = infraspecificEpithet)

#Giving an unique ID number for each record
paep_gbif <- cbind(id = 1:nrow(paep_gbif), paep_gbif)

#Encoding strings as UTF-8 (R version 4)
#class(paep_gbif) <- "data.frame" #if the current class is "data.table" "data.frame",
#the following code does not work properly. 
#for(i in 1:ncol(paep_gbif)){
#if(is.character(paep_gbif[ , i])){
#   Encoding(paep_gbif[ , i]) <- "UTF-8" 
#  }
#}

#=============#
# speciesLink #
#=============#

#Reading spLink 
paep_spLink <- fread(file = "datasets/Paepalanthus/speciesLink_paepalanthus_33895_20200826174927.txt", 
                     na.strings = c("", NA), stringsAsFactors = FALSE, encoding = "UTF-8")

#Selecting important attributes
paep_spLink <- paep_spLink %>% dplyr::select(institutioncode,
                                      collectioncode,
                                      catalognumber,
                                      genus,
                                      species,
                                      subspecies, 
                                      basisofrecord,
                                      identifiedby,
                                      collector,
                                      collectornumber,
                                      stateprovince,
                                      county,
                                      locality,
                                      longitude,
                                      latitude)

#Coercing coords into numeric values (NA's are introduced by coercion in observations
#with coordinate values as 'Bloqueada')
paep_spLink$longitude <- as.numeric(as.character(paep_spLink$longitude))
paep_spLink$latitude <- as.numeric(as.character(paep_spLink$latitude))

#Giving an unique ID number for each record
paep_spLink <- cbind(id = (nrow(paep_gbif) + 1):(nrow(paep_gbif) + nrow(paep_spLink)), paep_spLink)

#Encoding strings as UTF-8 
#class(paep_spLink) <- "data.frame" #if the current class is "data.table" "data.frame",
#the following code does not work properly. 
#for(i in 1:ncol(paep_spLink)){
#if(is.character(paep_spLink[ , i])){
# Encoding(paep_spLink[ , i]) <- "UTF-8" 
#}
#}

#==============================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (35,393), and adding a column to define the original dataset
#for each observation. 
merge.with.source <- function(x, y, name.x = "X", name.y = "Y") {
  x.df <- cbind(x, datsrc.x = name.x)
  y.df <- cbind(y, datsrc.y = name.y)
  merged.df <- merge(x = x.df,
                     y = y.df,
                     all = TRUE)
  merged.df[is.na(merged.df$datsrc.x), "datsrc.x"] <- ""
  merged.df[is.na(merged.df$datsrc.y), "datsrc.y"] <- ""
  merged.df$datsrc <- paste(merged.df$datsrc.x, merged.df$datsrc.y, sep = "")
  merged.df$datsrc.x <- rm()
  merged.df$datsrc.y <- rm()
  return(merged.df)
}

paep <- merge.with.source(x = paep_gbif,
                          y = paep_spLink,
                          name.x = "gbif",
                          name.y = "splink")

rm(merge.with.source, paep_gbif, paep_spLink)

#==============================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (observations with NA for any of those attributes
#were not considered) (28,740)
a <- paep
a.prime <- a[!is.na(a$institutioncode) &
               !is.na(a$collectioncode) &
               !is.na(a$catalognumber), ]
a.na <- a[is.na(a$institutioncode) |
            is.na(a$collectioncode) |
            is.na(a$catalognumber), ]
a <- unique(a.prime, by = c("institutioncode", 
                            "collectioncode",
                            "catalognumber"))
a <- rbind(a, a.na)
paep <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
paep <- paep %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
paep$latitude[paep$latitude == 0] <- NA
paep$longitude[paep$longitude == 0] <- NA

#Removing registers without identifier name (15,310)
plyr::count(paep$identifiedby)
paep$identifiedby[paep$identifiedby %in% c("?", "-", "#", "##", "//05/1987",
                                           "24/X/2013", "2015") | paep$identifiedby == 0] <- NA
paep <- paep %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- paep$stateprovince
#lookup_states <- data.frame("Incorreto" = plyr::count(province)[ , 1])
lookup_states <- fread(file = "lists/lookup_states.csv", na.strings = c("", NA), encoding = "UTF-8")
#write.csv(lookup_states, file = "lists/lookup_states.csv", row.names = F)
get_states <- lookup_states$Incorreto
names(get_states) <- lookup_states$Correto
for(i in 1:length(province)){
  for(j in 1:length(get_states)){
    if(province[i] == unname(get_states[j]) & !is.na(province[i])){
      province[i] <- names(get_states[j])
    }
  }
}
paep$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
paep$stateprovince[paep$stateprovince == "?" | paep$stateprovince == "-"] <- NA
plyr::count(paep$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur 
#according to Silveira et al (2016) (11,738). NA's are not considered.
paep <- paep %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                    "Minas Gerais", 
                                                                    "Bahia", 
                                                                    "Pernambuco", 
                                                                    "Paraiba", 
                                                                    "Mato Grosso",
                                                                    "Distrito Federal"))

#Removing records without species level identification (11,284)
paep <- paep %>% filter(!is.na(species))
paep <- paep %>% filter(!species %in% c("sp.", "sp1"))
#plyr::count(paep$species)

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (11,257)
#plyr::count(paep$basisofrecord)
paep <- paep %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                             "s", "S", "PreservedSpecimen"))
paep <- paep %>% select(-basisofrecord)

#==============================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

#Creating a vector to work with
identifiedby <- paep$identifiedby

#Generating a list of identifiers ordered by frequency of identifications (more identifications
#is supposed to be related with a better identification)
#Counting check
#names.count <- as.data.frame(plyr::count(identifiedby))
#names.count <- names.count[order(-names.count$freq), ]

#H. N. Moldenke
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("H. N.", "Moldenke"),
                              replace.by = "H. N. Moldenke",
                              not.replace = c(),
                              replace = c("Meldenke","Harold N. Moldenke",
                                          "H.L. Mello Barreto; Moldenke"))

#N. Hensold
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("N.", "Hensold"),
                              replace.by = "N. Hensold",
                              not.replace = c(),
                              replace = c("Nancy Hensold","N. Hensold : Field Museum of Natural History"))

#A. M. Giulietti
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("A. M.", "Giulietti"),
                              replace.by = "A. M. Giulietti",
                              not.replace = c(),
                              replace = c("Sano, P.T. & Giulietti, ªM.",
                                          "A.M. Giulietti; L.R. Parra","A.M. Giulietti & L.R. Parra",
                                          "Giuletti A.M. & Parra L.R.","A.M.Giulietti & M.J.G.Andrade",
                                          "Ana Maria Giuliti","A.Giulietti; L.P. Lazzari",
                                          "Giulietti, AM; Miranda Silva, EB","A.M.Giulietti & E.B.Miranda",
                                          "A.M. Giulietti & N. Hensold","A.M. Giulietti; L.P. Lazzari",
                                          "Guilietti","M. I. G. Andrade, A. M. Giulietti" ))

#M. L. O. Trovó
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("M. L. O.", "Trovó"),
                              replace.by = "M. L. O. Trovó",
                              not.replace = c(),
                              replace = c("M.Trovó","Trovó, MLO; Sano, PT; Echternacht, L",
                                          "M.I.O. Trovão", "M.C.O.Trovo","M. Trovó & C. Sarquis",
                                          "Trovó, MLO; Echternacht, L","Trovo#?#, M.",
                                          "M. Trovo","Trovo","Trauo, M.I.O.","M.I.O. Trouó",
                                          "M. Trovó & R. Ramos","M. Trovó 2016","M. Trovó & C.O. Andrino" ,
                                          "M.L.O.Trovó & C.Sarquis","M.L.O. Trovó & C. Sarquis",
                                          "M.L.O. Trovó & C. Sarquis","M.L.O. Trovó & P.T. Sano",
                                          "M.L.O. Trovó/04-XI-2014", "Trovó, MLO; Echternacht, LA",
                                          "M.L.O. Trovó; L. Sauthier","M.L.O. Trovó; R. Ramos",
                                          "Trovo & Sano","Marcelo Trovó","M. Trovo (SPF) 2008-10-22",
                                          "M.L.O. Trovó; C. Sarquis","Sarquis, C.; Trovó, M.",
                                          "Marcelo Trovó Lopes de Oliveira"))

#L. Echternacht
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("L.", "Echternacht"),
                              replace.by = "L. Echternacht",
                              not.replace = c(),
                              replace = c("A. Diaz; L. Echternacht",
                                          "Lívia Echternacht","L.Echternacht & M.Trovó",
                                          "Livia Echternacht","L.E. Andrade","Andrade, E.L."))

#P. T. Sano
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("P. T.", "Sano"),
                              replace.by = "P. T. Sano",
                              not.replace = c("Santos, R.M.","P. C. Standley",
                                              "D. T. Lano"),
                              replace = c("Sano, P.T.; Costa, F.N.; Trovó, M.L.O.; Echternacht, L.; Andrino, C.",
                                          "P.T. Sano; L.R. Parra","Sano, PT; Trovó, MLO",
                                          "P.T. Sano V.","Takeo, P.","Paulo T. Sano"))

#C. O. Andrino
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("C. O.", "Andrino"),
                              replace.by = "C. O. Andrino",
                              not.replace = c(),
                              replace = c("C.Andrino & F.N.Costa"))

#F. N. Costa
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("F. N.", "Costa"),
                              replace.by = "F. N. Costa",
                              not.replace = c("Costa, C.J.","Costa e Silva, M.B."),
                              replace = c("Costa, FN; Andrino, CO","F.N. Costa & . Andrino",
                                          "F.N.Costa & C.Adrino","F.N.Costa & C.Andrino",
                                          "F.N.Costa & C.O.Andrino"))
#B. Mourão
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("B.", "Mourão"),
                              replace.by = "B. Mourão",
                              not.replace = c(),
                              replace = c())

#M. J. G. Andrade
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M. J. G.", "Andrade"),
                              replace.by = "M. J. G. Andrade",
                              not.replace = c("Andrade, T.","Andrade, M.M.",
                                              "Andrade, MC","Andrade M.C."),
                              replace = c("M.J.G.Andrade & A.M.Giulietti"))

#E. B. Miranda
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("E. B.", "Miranda"),
                              replace.by = "E. B. Miranda",
                              not.replace = c("Miranda, AM","Miranda, A.M." ),
                              replace = c("E.Miranda-Silva"))

#L. J. Sauthier
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("L. J.", "Sauthier"),
                              replace.by = "L. J. Sauthier",
                              not.replace = c(),
                              replace = c())

#S. Splett
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("S.", "Splett"),
                              replace.by = "S. Splett",
                              not.replace = c(),
                              replace = c())

#M. M. Unwin
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M. M.", "Unwin"),
                              replace.by = "M. M. Unwin",
                              not.replace = c(),
                              replace = c())

#V. L. Scatena
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("V. L.", "Scatena"),
                              replace.by = "V. L. Scatena",
                              not.replace = c(),
                              replace = c())

#A. L. Silva
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("A. L.", "Silva"),
                              replace.by = "A. L. Silva",
                              not.replace = c("Silva, R.R.","Silva, E.M.",
                                              "Silva, BG","Silva, EBM","Silveira, A.",
                                              "A. L. R. Oliveira","Silva, M.G. da",
                                              "Silva, J.","A. C. Servilha","A. Lima",
                                              "Silva, M.C."),
                              replace = c())

#M. T. C. Watanabe
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M. T. C.", "Watanabe"),
                              replace.by = "M. T. C. Watanabe",
                              not.replace = c(),
                              replace = c())

#H. M. L. Tissot-Squalli
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("H. M. L.", "Tissot-Squalli"),
                              replace.by = "H. M. L. Tissot-Squalli",
                              not.replace = c(),
                              replace = c("M.L.Tissot-Squalli H."))

#R. Ramos
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("R.", "Ramos"),
                              replace.by = "R. Ramos",
                              not.replace = c("Ramos, A.E."),
                              replace = c("R.Ramos & R.B.Almeida"))

#L. R. Parra
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("L. R.", "Parra"),
                              replace.by = "L. R. Parra",
                              not.replace = c(),
                              replace = c("Pana, L.R. & Sano, P.T." ))

#W. Ruhland
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("W.", "Ruhland"),
                              replace.by = "W. Ruhland",
                              not.replace = c(),
                              replace = c())

#F. Koernicke
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("F.", "Koernicke"),
                              replace.by = "F. Koernicke",
                              not.replace = c(),
                              replace = c("Kornicke"))
#L. B. Smith
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("L. B.", "Smith"),
                              replace.by = "L. B. Smith",
                              not.replace = c("S. F. Smith"),
                              replace = c("Lyman B. Smith.","Lyman B. Smith"))

#Replacing column
paep$identifiedby <- identifiedby

#Filtering (9,578)
specialists <- c("H. N. Moldenke",
                 "N. Hensold",
                 "A. M. Giulietti",
                 "M. L. O. Trovó",
                 "L. Echternacht",
                 "P. T. Sano",
                 "C. O. Andrino",
                 "F. N. Costa",
                 "B. Mourão",
                 "M. J. G. Andrade",
                 "E. B. Miranda",
                 "L. J. Sauthier",
                 "S. Splett",
                 "M. M. Unwin",
                 "V. L. Scatena",
                 "A. L. Silva",
                 "M. T. C. Watanabe",
                 "H. M. L. Tissot-Squalli",
                 "R. Ramos",
                 "L. R. Parra",
                 "W. Ruhland",
                 "F. Koernicke",
                 "L. B. Smith")
paep <- paep %>% filter(identifiedby %in% specialists)

rm(identifiedby, specialists)

#====================================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

#Generating a column for scientific names (without authors and including infraspecific epithet)
paep$gen_sp <- paste(paep$genus,
                     paep$species,
                     paep$subspecies,
                     sep = " ")

#taxa <- plyr::count(paep$gen_sp)
#taxa$x <- as.character(taxa$x)
#taxa_suggested <- tibble("taxa" = taxa$x, "suggested" = NA)

#Removing NA (character derived from 'subspecies' attribute)
#for(i in 1:nrow(taxa_suggested)){
#taxa_suggested[i, 1] <- gsub("NA", "", taxa_suggested[i, 1])
#taxa_suggested[i, 1] <- trimws(taxa_suggested[i, 1])
#}

#Suggesting with flora
#for(i in 1:nrow(taxa_suggested)){
#taxa_suggested$suggested[i] <- suggest.names(taxa_suggested$taxa[i])
#}

#Writing *.csv for manual checking
#write.csv(taxa_suggested, file = "lists/Paepalanthus/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Paepalanthus/taxa_corrected.csv", stringsAsFactors = F, 
                           na.strings = c("NA",""))

#lookup 
taxa_gensp <- tibble(gen = NA, sp = NA, var = NA, .rows = nrow(taxa_corrected))
for(i in 1:nrow(taxa_corrected)){
  str <- strsplit(taxa_corrected$corrected[i], split = " ")[[1]]
  if(length(str) == 2){
    taxa_gensp$gen[i] <- str[1]
    taxa_gensp$sp[i] <- str[2]
  } else if(length(str) == 4){
    taxa_gensp$gen[i] <- str[1]
    taxa_gensp$sp[i] <- str[2]
    taxa_gensp$var[i] <- paste(str[3], str[4], sep = " ")
  }
}
taxa_gensp$replace <- taxa_corrected$taxa

#Removing invalid taxa (9,545)
paep$gen_sp <- gsub("NA", "", paep$gen_sp)
paep$gen_sp <- trimws(paep$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen)]
paep <- paep[!paep$gen_sp %in% invalid_taxa, ]

#Correcting the dataset (maintaining vars)
taxa_gensp$replace <- taxa_corrected$taxa
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(paep)){
  for(j in 1:nrow(taxa_gensp)){
    if(paep$gen_sp[i] == taxa_gensp$replace[j]){
      paep$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], taxa_gensp$var[j], sep = " ")
    }
  }
}

#Checking if there is any last corrections to be made
paep$gen_sp[!paep$gen_sp %in% paste(taxa_gensp$gen, taxa_gensp$sp, taxa_gensp$var,sep = " ")]

#Replacing NA by ""
paep$gen_sp <- gsub("NA", "", paep$gen_sp)
paep$gen_sp <- trimws(paep$gen_sp)

rm(taxa_corrected, taxa_gensp, invalid_taxa, str)

#======================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (4,448)
paep_coord <- paep %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (5,097)
paep_noCoord <- paep %>% filter(is.na(latitude) | is.na(longitude))

#======================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (4,262)
paep_coordFlagged <- paep_coord %>% clean_coordinates(lon = "longitude",
                                                      lat = "latitude",
                                                      species = "gen_sp",
                                                      value = "flagged",
                                                      tests = c("equal", "gbif", 
                                                                "institutions", 
                                                                "outliers", "seas",
                                                                "zeros"))

invalid_coords <- paep_coord[paep_coordFlagged == FALSE, ]
paep_coordClean <- paep_coord[paep_coordFlagged  == TRUE, ]

#Binding invalid_coords to paep_noCoord (8,146)
paep_noCoord <- rbind(paep_noCoord, invalid_coords)

rm(invalid_coords, paep_coordFlagged, paep_coord)

#Writing *.csv 
#write.csv(paep_coordClean, file = "Datasets/Paepalanthus/paep_coordClean.csv", row.names = FALSE)
#======================================================================================#

#===============================#
# CLEANING MUNICIPALITIES NAMES #
#===============================#

#Dividing dataset in registers with and without values for stateprovince (1,130 and 7,016)
paep_noState <- paep_noCoord %>% filter(is.na(stateprovince))
paep_state <- paep_noCoord %>% filter(!is.na(stateprovince))

#NA values for municipalities (2,380)
paep_na <- paep_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(paep_noState$county))
#names.count <- names.count[order(-names.count$freq), ]

#Reading list of municipalities according to CTB (Divisão Territorial Brasileira - IBGE) - ftp://geoftp.ibge.gov.br/organizacao_do_territorio/estrutura_territorial/divisao_territorial/
mun_concla <- fread(file = "lists/RELATORIO_DTB_BRASIL_MUNICIPIO.csv",
                    na.strings = c("", NA), stringsAsFactors = FALSE, encoding = "UTF-8")
mun_concla <- mun_concla %>% filter(Nome_UF %in% c("Bahia", "Goiás",
                                                   "Mato Grosso",
                                                   "Minas Gerais",
                                                   "Pernambuco",
                                                   "Paraíba",
                                                   "Distrito Federal"))

#Defining projection
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

#Reading shapefiles
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
#br <- readOGR("shapefiles/br_unidades_da_federacao/BRUFE250GC_SIR.shp") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm
mun <- readOGR("shapefiles/br_municipios/BRMUE250GC_SIR.shp") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#Projecting
proj4string(cr) <- crswgs84
#br <- spTransform(br, crswgs84)
mun <- spTransform(mun, crswgs84)

#Defining municipalities in which the campos rupestres occur
list_mun <- over(cr, mun)
list_mun <- unique(list_mun$NM_MUNICIP)
list_mun <- as.character(list_mun)
#Encoding(list_mun) <- "UTF-8" #Run this line if R fails to encode the vector created above

#Testing intersection
#plot(mun[mun$NM_MUNICIP %in% list_mun, ]) 

#Which municipalities are homonyms?
homonyms <- mun$NM_MUNICIP[mun$NM_MUNICIP %in% list_mun]
homonyms <- correct.mun(homonyms[duplicated(homonyms)])

#Standardizing municipalities names - lists and data sets
list_mun_std <- correct.mun(list_mun)
mun_concla$Nome_Município <- correct.mun(mun_concla$Nome_Município)
paep_noState$municipality_gbif_std <- correct.mun(paep_noState$municipality_gbif)
paep_state$municipality_gbif_std <- correct.mun(paep_state$municipality_gbif)
paep_noState$county_std <- correct.mun(paep_noState$county)
paep_state$county_std <- correct.mun(paep_state$county)

#Selecting mun_concla municipalities in which the campos rupestres occur
mun_concla_cr <- mun_concla[mun_concla$Nome_Município %in% list_mun_std, ]

#nrow(mun_concla_cr) is not equal to length(list_mun_std). Why?
#Municipalities from the São Paulo state
#list_mun_std[!list_mun_std %in% mun_concla_cr$Nome_Município]

#Manually checking and correcting by crossing the list assigned as 'check' with the mun_concla_cr dataset. 
#check <- data.frame(plyr::count(paep_noState$municipality_gbif_std)) #ok
paep_noState$municipality_gbif_std[paep_noState$municipality_gbif_std == 
                                     "Anaje"] <- "Anage"
paep_noState$municipality_gbif_std[paep_noState$municipality_gbif_std == 
                                     "Gouvea"] <- "Gouveia"
paep_noState$municipality_gbif_std[paep_noState$municipality_gbif_std == 
                                     "Sao Joao D'Alianca"] <- "Sao Joao D'alianca"
paep_noState$municipality_gbif_std[paep_noState$municipality_gbif_std 
                                   %in% c("Brasil", "Goias",
                                          "Mun?")] <- NA

#check <- data.frame(plyr::count(paep_noState$county_std)) #ok
paep_noState$county_std[paep_noState$county_std 
                        %in% c("Brasil", "5Km Sul da Cidade")] <- NA
paep_na <- rbind(paep_na, 
                 paep_noState[which(is.na(paep_noState$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                      is.na(paep_noState$county_std)), 1:18]) 
paep_noState <- paep_noState[-which(is.na(paep_noState$municipality_gbif_std) & #removing registers with NA values for these two attributes from paep_noState
                                      is.na(paep_noState$county_std)), ]

#check <- data.frame(plyr::count(paep_state$municipality_gbif_std))
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std %in% 
                                   c("Alto Garca")] <- "Alto Garcas"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std %in% 
                                   c("Alto Paraiso","Alto Paraiso Goias")] <- "Alto Paraiso de Goias"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Abaira   Piata"] <- "Abaira"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Anaje"] <- "Anage"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Barra de Mendes"] <- "Barra do Mendes"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Betania/Floresta"] <- "Betania"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std %in% 
                                   c("Braslilia")] <- "Brasilia"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Catugi"] <- "Catuji"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Cajazeiras Sao Jose das Piranhas"] <- "Cajazeiras"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Cidade Ecletica"] <- "Santo Antonio do Descoberto"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std %in% 
                                   c("Conceicao do Mato de Dentro",
                                     "Conc do Mato Dentro")] <- "Conceicao do Mato Dentro"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Corumba  de Goias"] <- "Corumba de Goias"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std %in% 
                                   c("Delfinopolis (?)","Delfinopoilis")] <- "Delfinopolis"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Felixandia"] <- "Felixlandia"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Gouvea"] <- "Gouveia"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std %in% 
                                   c("Jaboticatuba")] <- "Jaboticatubas"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Olhos D#?#Agua"] <- "Olhos D'agua"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Planaltina de Goias"] <- "Planaltina"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std 
                                 %in% c("Presidente Kubitchek")] <- "Presidente Kubitschek"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Rio das Contas"] <- "Rio de Contas"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Santana do Pirapama"] <- "Santana de Pirapama"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Sao Tome das Letras"] <- "Sao Thome das Letras"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Sao Joao da Alianca"] <- "Sao Joao D'alianca"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std 
                                 %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                        "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std == 
                                   "Terezina de Goias"] <- "Teresina de Goias"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std 
                                 %in% c("Varzea de Palma")] <- "Varzea da Palma"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std 
                                 %in% c("Vila Bela Santissima Trindade", "Vila Bela da Sma Trindade",
                                        "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std 
                                 %in% c("Chapada da Cotagem", "Chapada Diamantina",
                                        "Chapada do Araripe", "Chapada dos Guimaraes",
                                        "Chapada dos Veadeiros", "Paraiba", "Pernambuco",
                                        "Serra do Cipo", "Serra do Espinhaco", "Serra do Tombador",
                                        "Ufba Ondina", "Veadeiros", "Bahia", "Chapadao do Ceu", "Goias", 
                                        "Minas Gerais", "Serra do Cipo", "Serra do Espinhaco", 
                                        "Serra do Tombador", "Veadeiros", "Bahia", "Betim/Brumadinho",
                                        "Coletada A Margem da Estrada de Macambinha","Br 135 Km 404",
                                        "Lavras Sao Joao Del Rey", "Ouro Preto/ Mariana")] <- NA


#check <- data.frame(plyr::count(paep_state$county_std))
paep_state$county_std[paep_state$county_std %in% 
                        c("Alto Paraiso", "Alto Paraa-so de Goias")] <- "Alto Paraiso de Goias"
paep_state$county_std[paep_state$county_std == 
                        "Anaje"] <- "Anage"
paep_state$county_std[paep_state$county_std == 
                        "Barra de Mendes"] <- "Barra do Mendes"
paep_state$county_std[paep_state$county_std == 
                        "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std %in% 
                                   c("Conceicao do Mato de Dentro",
                                     "Conc do Mato Dentro","Conceicao Doato Dentro")] <- "Conceicao do Mato Dentro"
paep_state$county_std[paep_state$county_std == 
                        "Cristalina Mun"] <- "Cristalina"
paep_state$county_std[paep_state$county_std %in% 
                        c("Golvea")] <- "Golveia"
paep_state$county_std[paep_state$county_std %in% 
                        c("Grao Mogol Mun")] <- "Grao Mogol"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std %in% 
                                   c("Jaboticatuba")] <- "Jaboticatubas"
paep_state$county_std[paep_state$county_std == 
                        "Joaquim Fela-cio"] <- "Joaquim Felicio"
paep_state$county_std[paep_state$county_std %in% 
                        c("Morro do Chapeu Mun")] <- "Morro do Chapeu"
paep_state$county_std[paep_state$county_std == 
                        "Rio das Contas"] <- "Rio de Contas"
paep_state$county_std[paep_state$county_std %in% 
                        c("Sao Joao da Alianca", 
                          "Sao Joao D##Alianca")] <- "Sao Joao D'alianca"
paep_state$county_std[paep_state$county_std %in% 
                        c("Sao Tome das Letras ")] <- "Sao Thome das Letras "
paep_state$county_std[paep_state$county_std %in% 
                        c("Vila Bela de Santa-ssima Trindade", 
                          "Vila Bela de Santissima Trindade",
                          "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
paep_state$county_std[paep_state$county_std == 
                        "Santana de Pirapama"] <- "Santana do Pirapama"
paep_state$municipality_gbif_std[paep_state$municipality_gbif_std 
                                 %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                        "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
paep_state$county_std[paep_state$county_std == 
                        "Brasa-lia"] <- "Brasilia"
paep_state$county_std[paep_state$county_std 
                      %in% c("", "Chapada dos Guimaraes",
                             "Chapada dos Veadeiros", "Chapada Gaucha",
                             "Goias", "Minas Gerais", "No Disponible",
                             "Pe", "Ufba Ondina", "Veadeiros", 
                             "Minas Gerais", "Mato Grosso", "Serra de Ibitipoca",
                             "Serra do Cabral", "Serra do Cipo", "Serra do Espinhaco",
                             "Entre Serro E Lagoa Santa","Lavras Sao Joao Del Rey")] <- NA
paep_na <- rbind(paep_na, 
                 paep_state[which(is.na(paep_state$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                    is.na(paep_state$county_std)), 1:18]) 
paep_state <- paep_state[-which(is.na(paep_state$municipality_gbif_std) & #removing registers with NA values for these two attributes from paep_noState
                                  is.na(paep_state$county_std)), ]

#Filtering paep_noState and paep_state with mun_concla_cr
paep_noState_filt <- paep_noState %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                                county_std %in% mun_concla_cr$Nome_Município)
paep_state_filt <- paep_state %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                            county_std %in% mun_concla_cr$Nome_Município)

#===================================================================================================#

#================================#
# INFERRING MUNICIPALITIES NAMES #
#================================#

#Standardizing 'locality'
paep_na$locality_std <- correct.mun(paep_na$locality)

#Vector with municipalities names from mun_concla to be used with grepl
grepl_munc_concla_cr <- c()
for(i in 1:nrow(mun_concla_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", mun_concla_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(paep_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(paep_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = paep_na$locality_std[j])){
      vec[j] <- mun_concla_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(paep_na$locality_std)] <- NA

#New column with inferred municipality names
paep_na$municipality <- vec

#Removing registers with NA for 'municipality'
paep_na <- paep_na %>% filter(!is.na(municipality))

#Concatenating information on municipality into a unique column 
paep_noState_filt$municipality <- NA
for(i in 1:nrow(paep_noState)){
  if(!is.na(paep_noState_filt$municipality_gbif_std[i])){
    paep_noState_filt$municipality[i] <- paep_noState_filt$municipality_gbif_std[i] 
  } else if(!is.na(paep_noState_filt$county_std[i])){
    paep_noState_filt$municipality[i] <- paep_noState_filt$county_std[i]
  }
}

paep_state_filt$municipality <- NA
for(i in 1:nrow(paep_state)){
  if(!is.na(paep_state_filt$municipality_gbif_std[i])){
    paep_state_filt$municipality[i] <- paep_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(paep_state_filt$county_std[i])){
    paep_state_filt$municipality[i] <- paep_state_filt$county_std[i]
  }
}

#Concatenating data sets
paep_noState_filt <- paep_noState_filt[ , -which(names(paep_noState_filt) 
                                                 %in% c("county","county_std",
                                                        "municipality_gbif",
                                                        "municipality_gbif_std"))]

paep_state_filt <- paep_state_filt[ , -which(names(paep_state_filt) 
                                             %in% c("county","county_std",
                                                    "municipality_gbif",
                                                    "municipality_gbif_std"))]
paep_na <- paep_na[ , -which(names(paep_na) 
                             %in% c("county", "municipality_gbif", "locality_std"))]

paep_noCoord_inf<- rbind(paep_noState_filt, paep_state_filt, paep_na)

#Registers occurring in homonyms municipalities
reg_hom <- paep_noCoord_inf[paep_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
paep_noCoord_inf <- paep_noCoord_inf %>% filter(!municipality %in% homonyms)
paep_noCoord_inf <- rbind(paep_noCoord_inf, reg_hom) #binding them after removing those without state/province information

rm(paep_na, paep_noState, paep_noState_filt, paep_state, paep_state_filt,
   mun_concla, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, paep_noCoord, reg_hom)

#=====================================================================================================#

#===============================================#
# INFERRING COORDINATES BASED ON 'municipality' #
#===============================================#

#Loading centroids' data set
mun_centr <- fread(file = "lists/Neotropics_municipalities.csv", na.strings = c("", NA))

#Filtering for Brazilian municipalities
mun_centr_br <- mun_centr %>% filter(country == "Brazil")

#Standardizing municipalities names
mun_centr_br$`minor area` <- correct.mun(mun_centr_br$`minor area`)

#Filtering for cr municipalities
mun_centr_cr <- mun_centr_br %>% filter(`minor area` %in% mun_concla_cr$Nome_Município &
                                          `major area` %in% mun_concla_cr$Nome_UF)

#Removing homonyms from the data set and replacing them 
centr_hom <- mun_centr_br %>% filter(`minor area` %in% homonyms)
mun_centr_cr <- mun_centr_cr %>% filter(!`minor area` %in% homonyms)

#Some municipalities are not included in the data set for centroids
#not_included <- mun_concla_cr %>% filter(!Nome_Município %in% mun_centr_cr$`minor area`)

#Inferring coordinates based on 'municipality'
for(i in 1:length(paep_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(paep_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      paep_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      paep_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Homonyms
for(i in 1:length(paep_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(paep_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       paep_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(paep_noCoord_inf$stateprovince[i])){
      paep_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      paep_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing NA's for coordinates (3,239)
paep_noCoord_inf <- paep_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, mun_concla_cr, centr_hom, homonyms)

#=================================================================================================#

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unnimportant columns (7,501)
paep_noCoord_inf <- paep_noCoord_inf %>% select(id,
                                                institutioncode,
                                                collectioncode,
                                                catalognumber,
                                                gen_sp,
                                                subspecies,
                                                identifiedby,
                                                latitude,
                                                longitude)

paep_coordClean <- paep_coordClean %>% select(id,
                                              institutioncode,
                                              collectioncode,
                                              catalognumber,
                                              gen_sp,
                                              subspecies,
                                              identifiedby,
                                              latitude,
                                              longitude)

paep_clean <- rbind(paep_coordClean, paep_noCoord_inf)

#Is there any duplicated registers? 
#paep_clean[duplicated(paep_clean$id) == TRUE, ] #no

rm(paep, paep_coordClean, paep_noCoord_inf, i, j, correct.mun, generate.names,
   replace.names, titling)

#Writing *.csv 
#write.csv(paep_clean, file = "Datasets/Paepalanthus/paep_clean.csv", row.names = FALSE)

#====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Reading clean dataset
paep_clean <- read.csv(file = "datasets/Paepalanthus/paep_clean.csv", na.strings = c("", NA))

#Standardizing gen_sp column
paep_clean$gen_sp <- gsub(" ", "_", paep_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the cr shapefile and, then, with the grids 
coords <- paep_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- paep_clean
coords$id_2 <- 1:nrow(coords)
coords_2 <- coords_2 %>% filter(!is.na(ID))
coords <- coords %>% filter(id_2 %in% coords_2$id_2)
coords_2 <- coords 
coordinates(coords_2) <- ~ longitude + latitude
proj4string(coords_2) <- crswgs84
coords_2 <- over(coords_2, grids_cr)
coords_2$id_2 <- 1:nrow(coords_2) 
coords$id_2 <- 1:nrow(coords_2)
coords_2 <- coords_2 %>% filter(!is.na(id))
coords <- coords %>% filter(id_2 %in% coords_2$id_2)
coords$id_grid <- coords_2$id
coords <- coords[ , -which(colnames(coords) == "id_2")]

#Saving dataset
write.csv(coords, "datasets/Paepalanthus/paep_cr.csv", row.names = FALSE)
