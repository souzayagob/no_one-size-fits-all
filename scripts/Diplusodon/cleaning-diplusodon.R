#This code is intended to clean datasets downloaded from GBIF and speciesLink, so they contain
#only records from the campos rupestres

#====================================================================================================#

library(tidyverse) 
library(data.table)
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 
library(raster)

source('~/Dropbox/Mestrado/Bancos_organizing/Scripts/functions.R')

#====================================================================================================#

#==================#
# READING DATASETS #
#==================#

#======#
# GBIF #
#======#

#Reading gbif data
diplusodon_gbif <- fread(file = "Datasets/Diplusodon/0008173-190415153152247_diplusodon_gbif/occurrence.txt",
                     na.strings = c("", NA))

#Selecting important attributes
diplusodon_gbif <- diplusodon_gbif %>% dplyr::select(institutionCode,
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
diplusodon_gbif <- diplusodon_gbif %>% rename("species" = specificEpithet,
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
diplusodon_gbif <- cbind(id = 1:nrow(diplusodon_gbif), diplusodon_gbif)

#=============#
# speciesLink #
#=============#

#Reading spLink data
diplusodon_spLink <- fread(file = "Datasets/Diplusodon/diplusodon_spLink.txt", 
                       na.strings = c("", NA))

#Selecting important attributes
diplusodon_spLink <- diplusodon_spLink %>% dplyr::select(institutioncode,
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
#with coordinate values equal to 'Bloqueada')
#diplusodon_spLink$longitude <- as.numeric(as.character(diplusodon_spLink$longitude))
#diplusodon_spLink$latitude <- as.numeric(as.character(diplusodon_spLink$latitude))

#Giving an unique ID number for each record
diplusodon_spLink <- cbind(id = (nrow(diplusodon_gbif) + 1):(nrow(diplusodon_gbif) + nrow(diplusodon_spLink)), diplusodon_spLink)

#=====================================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (18,497), and adding a column to define the original dataset
#for each observation
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

diplusodon <- merge.with.source(x = diplusodon_gbif,
                            y = diplusodon_spLink,
                            name.x = "gbif",
                            name.y = "splink")

rm(merge.with.source, diplusodon_gbif, diplusodon_spLink)

#=====================================================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (NA's not considered) (13,629)
a <- diplusodon
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
diplusodon <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
diplusodon <- diplusodon %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA (coords)
diplusodon$latitude[diplusodon$latitude == 0] <- NA
diplusodon$longitude[diplusodon$longitude == 0] <- NA

#Removing registers without identifier's name (8,301)
#plyr::count(diplusodon$identifiedby)
diplusodon$identifiedby[diplusodon$identifiedby %in% c("?", "20-5-")] <- NA
diplusodon <- diplusodon %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- diplusodon$stateprovince
#Writing lookup table
#write.csv(lookup_states, file = "lists/lookup_states.csv", row.names = F)
#lookup_states <- data.frame("Incorreto" = plyr::count(province)[ , 1])
#Reading lookup table
lookup_states <- fread(file = "lists/lookup_states.csv", na.strings = c("", NA))
get_states <- lookup_states$Incorreto
names(get_states) <- lookup_states$Correto
for(i in 1:length(province)){
  for(j in 1:length(get_states)){
    if(province[i] == unname(get_states[j]) & !is.na(province[i])){
      province[i] <- names(get_states[j])
    }
  }
}
diplusodon$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
#plyr::count(diplusodon$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur 
#according to Silveira et al (2016) (7,906)
diplusodon <- diplusodon %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                        "Minas Gerais", 
                                                                        "Bahia", 
                                                                        "Pernambuco", 
                                                                        "Paraiba", 
                                                                        "Mato Grosso",
                                                                        "Distrito Federal"))

#Removing records without species level identification (7,664)
diplusodon <- diplusodon %>% filter(!is.na(species))
diplusodon <- diplusodon %>% filter(!species %in% c("sp."))
#plyr::count(diplusodon$species)

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (7,640)
#plyr::count(diplusodon$basisofrecord)
diplusodon <- diplusodon %>% filter(basisofrecord %in% c("PRESERVED_SPECIMEN", "S"))
diplusodon <- diplusodon %>% dplyr::select(-basisofrecord)

#=====================================================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

identifiedby <- diplusodon$identifiedby

#Generating a list of identifiers ordered by frequency of identifications (more identifications
#is supposed to be related with a better identification)
#Counting check
#names.count <- as.data.frame(plyr::count(identifiedby))
#names.count <- names.count[order(-names.count$freq), ]

#T. B. Cavalcanti
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("T. B.", "Cavalcanti"),
                              replace.by = "T. B. Cavalcanti",
                              not.replace = c(),
                              replace = c("TACIANA B. CAVALCANTI", "T.B.CAVALCANTI (SPF) V. LORENTIG",
                                          "T.B.CAVALCANTI (SPF) V. A. LOURTUG", "T.B.CAVALCANTE (SPP) V.H.LORUTEG",
                                          "Ticiana B. Cavalcanti", "Taciana Barbsa Cavalcanti", "Tacaiana B. Cavalcanti",
                                          "T.B", "T.B. Cavalcanti & S. Graham","Taciana B. Calvacanti", 
                                          "T. Barbosa Cavalcanti 1999-04-15", "T. Barbosa", "Tacina Cavalcanti",
                                          "Taciana B. Cavalcanti", "Taciana Barbosa Cavalcanti",
                                          ".FTaciana Barbosa cavalcanti",
                                          "Taciana aBarbosa Cavalcanti", "Taciana arbosa Cavalcanti",
                                          "Taciana Barbosa", "taciana Barbosa Cavalcanti",
                                          "Taciana BArbosa Cavalcanti", "Taciana Barbosa Cavalcanti et al",
                                          "Taciana Barbosa Cavarcanti", "Tacianana Barbosa Cavalcanti", 
                                          "Tacianan Barbosa Cavalcanti",
                                          "Barbosa Cavalcanti, T"))

#A. Lourteig
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("A.", "Lourteig"),
                              replace.by = "A. Lourteig",
                              not.replace = c("A. Lowtey"),
                              replace = c())

#S. A. Graham
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("S. A.", "Graham"),
                              replace.by = "S. A. Graham",
                              not.replace = c("S. A. Mori"),
                              replace = c("Shirley Graham", "S. Graham (MO)", 
                                          "S. Graham (MO), Aug.", "S.A. Graham (MO)",
                                          "Shirley A. Graham", "S.A. Graham, 2003 (MO)",
                                          "Shirley A Graham (MO)"))

#E. Koehne
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("E.", "Koehne"),
                              replace.by = "E. Koehne",
                              not.replace = c(),
                              replace = c("B. A. E. Koehne",
                                          "B. A. E. Koehne 1901-04-05",
                                          "B. A. E. Koehne 1989"))

#names.count <- as.data.frame(plyr::count(identifiedby))
#names.count <- names.count[order(-names.count$freq), ]

#Replacing column
diplusodon$identifiedby <- identifiedby

#Filtering (6,274)
specialists <- c("T. B. Cavalcanti",
                 "A. Lourteig",
                 "S. A. Graham",
                 "E. Koehne")
diplusodon <- diplusodon %>% filter(identifiedby %in% specialists)

rm(identifiedby, specialists)

#=====================================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

#Generating a column for scientific names (without authors and including infraspecific epithet)
diplusodon$gen_sp <- paste(diplusodon$genus,
                          diplusodon$species,
                          diplusodon$subspecies,
                          sep = " ")

#taxa <- plyr::count(diplusodon$gen_sp)
#taxa$x <- as.character(taxa$x)
#taxa_suggested <- tibble("taxa" = taxa$x, "suggested" = NA)

#Removing NA (character derived from 'subspecies' attribute)
#for(i in 1:nrow(taxa_suggested)){
#taxa_suggested[i, 1] <- gsub("NA", "", taxa_suggested[i, 1])
#taxa_suggested[i, 1] <- trimws(taxa_suggested[i, 1])
#}

#Suggesting with the package 'flora' 
#for(i in 1:nrow(taxa_suggested)){
 #taxa_suggested$suggested[i] <- suggest.names(taxa_suggested$taxa[i])
#}

#Writing *.csv for manual checking
#write.csv(taxa_suggested, file = "lists/Diplusodon/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Diplusodon/taxa_corrected.csv", stringsAsFactors = F, 
                           na.strings = c("NA",""))

#trimws for every column
for(i in 1:ncol(taxa_corrected)){
  taxa_corrected[ , i] <- trimws(taxa_corrected[ , i])
}

#Dataset with different columns for different epithets
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

#Removing invalid taxa (6,258)
diplusodon$gen_sp <- gsub("NA", "", diplusodon$gen_sp)
diplusodon$gen_sp <- trimws(diplusodon$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen)]
diplusodon <- diplusodon[!diplusodon$gen_sp %in% invalid_taxa, ]

#Correcting the dataset
diplusodon$gen_sp <- trimws(gsub("NA", "", diplusodon$gen_sp))
taxa_gensp$replace <- taxa_corrected$taxa
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(diplusodon)){
  for(j in 1:nrow(taxa_gensp)){
    if(diplusodon$gen_sp[i] == taxa_gensp$replace[j]){
      diplusodon$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

rm(taxa_corrected, taxa_gensp, str, invalid_taxa)

#=====================================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (3,322)
diplusodon_coord <- diplusodon %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (2,936)
diplusodon_noCoord <- diplusodon %>% filter(is.na(latitude) | is.na(longitude))

#Two registers, ids 14240 e 14277, have invalid coordinates
#Removing them and transfering to diplusodon_noCoord
diplusodon_coord <- diplusodon_coord %>% filter(!id %in% c(14240, 14277))
diplusodon_noCoord <- rbind(diplusodon_noCoord, diplusodon[diplusodon$id %in% c(14240, 14277), ])

rm(diplusodon)
#=====================================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning
diplusodon_coordFlagged <- diplusodon_coord %>% clean_coordinates(lon = "longitude",
                                                        lat = "latitude",
                                                        species = "gen_sp",
                                                        value = "flagged",
                                                        tests = c("equal", "gbif", 
                                                                  "institutions", 
                                                                  "outliers", "seas",
                                                                  "zeros"))


invalid_coords <- diplusodon_coord[diplusodon_coordFlagged == FALSE, ]
diplusodon_coordClean <- diplusodon_coord[diplusodon_coordFlagged  == TRUE, ]

#Binding invalid_coords to diplusodon_noCoord
diplusodon_noCoord <- rbind(diplusodon_noCoord, invalid_coords)

rm(invalid_coords, diplusodon_coordFlagged, diplusodon_coord)

#Writing *.csv
write.csv(diplusodon_coordClean, 
          file = "Datasets/Diplusodon/diplusodon_coordClean.csv", row.names = FALSE)
#=====================================================================================================#

#===============================#
# CLEANING MUNICIPALITIES' NAMES #
#===============================#

#Dividing dataset in registers with and without values for 'stateprovince' (370 and 2,657)
diplusodon_noState <- diplusodon_noCoord %>% filter(is.na(stateprovince))
diplusodon_state <- diplusodon_noCoord %>% filter(!is.na(stateprovince))

#Dataset with NA values for municipalities (924)
diplusodon_na <- diplusodon_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(diplusodon_noState$county))
#names.count <- names.count[order(-names.count$freq), ]

#Reading list of municipalities according to CTB (Divisão Territorial Brasileira - IBGE) - ftp://geoftp.ibge.gov.br/organizacao_do_territorio/estrutura_territorial/divisao_territorial/
mun_concla <- fread(file = "lists/RELATORIO_DTB_BRASIL_MUNICIPIO.csv",
                    na.strings = c("", NA))
mun_concla <- mun_concla %>% filter(Nome_UF %in% c("Bahia", "Goiás",
                                                   "Mato Grosso",
                                                   "Minas Gerais",
                                                   "Pernambuco",
                                                   "Paraíba",
                                                   "Distrito Federal"))

#Defining projection
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

#Reading shapefiles for campos rupestres and municipalities
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
mun <- readOGR("shapefiles/br_municipios/BRMUE250GC_SIR.shp") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#Projecting
proj4string(cr) <- crswgs84
mun <- spTransform(mun, crswgs84)

#Defining municipalities in which the campos rupestres occur
list_mun <- over(cr, mun)
list_mun <- unique(list_mun$NM_MUNICIP)
list_mun <- as.character(list_mun) 

#Testing intersection
#plot(mun[mun$NM_MUNICIP %in% list_mun, ]) 

#Which municipalities are homonyms?
homonyms <- mun$NM_MUNICIP[mun$NM_MUNICIP %in% list_mun]
homonyms <- correct.mun(homonyms[duplicated(homonyms)])

#Standardizing municipalities' names - lists and datasets
list_mun_std <- correct.mun(list_mun)
mun_concla$Nome_Município <- correct.mun(mun_concla$Nome_Município)
diplusodon_noState$municipality_gbif_std <- correct.mun(diplusodon_noState$municipality_gbif)
diplusodon_state$municipality_gbif_std <- correct.mun(diplusodon_state$municipality_gbif)
diplusodon_noState$county_std <- correct.mun(diplusodon_noState$county)
diplusodon_state$county_std <- correct.mun(diplusodon_state$county)

#Selecting mun_concla municipalities in which the campos rupestres occur
mun_concla_cr <- mun_concla[mun_concla$Nome_Município %in% list_mun_std, ]

#Manually checking and correcting 
#check <- data.frame(plyr::count(diplusodon_noState$municipality_gbif_std)) #ok
diplusodon_na <- rbind(diplusodon_na, 
                   diplusodon_noState[which(diplusodon_noState$municipality_gbif_std %in% 
                                          c("Brasil", "Chapada dos Guimaraes")), 1:18])
diplusodon_state$municipality_gbif_std[diplusodon_state$municipality_gbif_std == 
                                     "Sao Joao D'Alianca"] <- "Sao Joao D'alianca"

#check <- data.frame(plyr::count(diplusodon_noState$county_std)) #ok

#check <- data.frame(plyr::count(diplusodon_state$municipality_gbif_std)) #ok
diplusodon_state$municipality_gbif_std[diplusodon_state$municipality_gbif_std %in% 
                                    c("Alto Paraiso", "Alto do Paraiso de Goias")] <- "Alto Paraiso de Goias"
diplusodon_state$municipality_gbif_std[diplusodon_state$municipality_gbif_std == 
                                    "Caitite"] <- "Caetite"
diplusodon_state$municipality_gbif_std[diplusodon_state$municipality_gbif_std == 
                                         "Gouvea"] <- "Gouveia"
diplusodon_state$municipality_gbif_std[diplusodon_state$municipality_gbif_std == 
                                         "Piranopolis"] <- "Pirenopolis"
diplusodon_state$municipality_gbif_std[diplusodon_state$municipality_gbif_std == 
                                         "Santana de Riacho"] <- "Santana do Riacho"
diplusodon_state$municipality_gbif_std[diplusodon_state$municipality_gbif_std == 
                                         "Sao Joao Del Rey"] <- "Sao Joao Del Rei"
diplusodon_state$municipality_gbif_std[diplusodon_state$municipality_gbif_std == 
                                         "Sao Tome das Letras"] <- "Sao Thome das Letras"
diplusodon_na <- rbind(diplusodon_na, 
                   diplusodon_state[which(diplusodon_state$municipality_gbif_std %in% 
                                        c("Chapada dos Guimaraes", "Chapada dos Veadeiros", 
                                          "Goias", "Minas Gerais", "Serra do Espinhaco")), 1:18])

#check <- data.frame(plyr::count(diplusodon_state$county_std)) #ok
diplusodon_state$county_std[diplusodon_state$county_std == 
                         "Alto Paraa-so de Goias"] <- "Alto Paraiso de Goias"
diplusodon_state$county_std[diplusodon_state$county_std == 
                         "Brasa-lia"] <- "Brasilia"
diplusodon_state$county_std[diplusodon_state$county_std == 
                              "Joaquim Fela-cio"] <- "Joaquim Felicio"
diplusodon_state$county_std[diplusodon_state$county_std == 
                              "Sao Tome das Letras"] <- "Sao Thome das Letras"

#Filtering diplusodon_noState and diplusodon_state with mun_concla_cr (71 and 1,409)
diplusodon_noState_filt <- diplusodon_noState %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                                  county_std %in% mun_concla_cr$Nome_Município)
diplusodon_state_filt <- diplusodon_state %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                              county_std %in% mun_concla_cr$Nome_Município)

#=====================================================================================================#

#================================#
# INFERRING MUNICIPALITIES NAMES #
#================================#

#Standardizing 'locality'
diplusodon_na$locality_std <- correct.mun(diplusodon_na$locality)

#Vector with municipalities names from mun_concla to be used with 'grepl'
grepl_munc_concla_cr <- c()
for(i in 1:nrow(mun_concla_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", mun_concla_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(diplusodon_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(diplusodon_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = diplusodon_na$locality_std[j])){
      vec[j] <- mun_concla_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(diplusodon_na$locality_std)] <- NA

#New column with inferred municipality names
diplusodon_na$municipality <- vec

#Removing registers with NA for 'municipality' (386)
diplusodon_na <- diplusodon_na %>% filter(!is.na(municipality))

#=====================================================================================================#

#Concatenating information on municipality into a unique column 
diplusodon_noState_filt$municipality <- NA
for(i in 1:nrow(diplusodon_noState)){
  if(!is.na(diplusodon_noState_filt$municipality_gbif_std[i])){
    diplusodon_noState_filt$municipality[i] <- diplusodon_noState_filt$municipality_gbif_std[i] 
  } else if(!is.na(diplusodon_noState_filt$county_std[i])){
    diplusodon_noState_filt$municipality[i] <- diplusodon_noState_filt$county_std[i]
  }
}

diplusodon_state_filt$municipality <- NA
for(i in 1:nrow(diplusodon_state)){
  if(!is.na(diplusodon_state_filt$municipality_gbif_std[i])){
    diplusodon_state_filt$municipality[i] <- diplusodon_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(diplusodon_state_filt$county_std[i])){
    diplusodon_state_filt$municipality[i] <- diplusodon_state_filt$county_std[i]
  }
}

#Concatenating data sets (1,866)
diplusodon_noState_filt <- diplusodon_noState_filt[ , -which(names(diplusodon_noState_filt) 
                                                   %in% c("county","county_std",
                                                          "municipality_gbif",
                                                          "municipality_gbif_std"))]
diplusodon_state_filt <- diplusodon_state_filt[ , -which(names(diplusodon_state_filt) 
                                               %in% c("county","county_std",
                                                      "municipality_gbif",
                                                      "municipality_gbif_std"))]
diplusodon_na <- diplusodon_na[ , -which(names(diplusodon_na) 
                               %in% c("county", "municipality_gbif", "locality_std"))]
diplusodon_noCoord_inf<- rbind(diplusodon_noState_filt, diplusodon_state_filt, diplusodon_na)

#Registers occurring in homonyms municipalities (1,861)
reg_hom <- diplusodon_noCoord_inf[diplusodon_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
diplusodon_noCoord_inf <- diplusodon_noCoord_inf %>% filter(!municipality %in% homonyms)
diplusodon_noCoord_inf <- rbind(diplusodon_noCoord_inf, reg_hom) #binding them after removing those without state/province information

rm(diplusodon_na, diplusodon_noState, diplusodon_noState_filt, diplusodon_state, diplusodon_state_filt,
   mun_concla, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, diplusodon_noCoord, reg_hom)

#====================================================================================================#

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

#Inferring coordinates based on 'municipality'
for(i in 1:length(diplusodon_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(diplusodon_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      diplusodon_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      diplusodon_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Homonyms
for(i in 1:length(diplusodon_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(diplusodon_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       diplusodon_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(diplusodon_noCoord_inf$stateprovince[i])){
      diplusodon_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      diplusodon_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing NA's for coordinates (1,861)
diplusodon_noCoord_inf <- diplusodon_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, mun_concla_cr, centr_hom, homonyms)

#=====================================================================================================#

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unnimportant columns (5,060)
diplusodon_noCoord_inf <- diplusodon_noCoord_inf %>% dplyr::select(id,
                                                  institutioncode,
                                                  collectioncode,
                                                  catalognumber,
                                                  gen_sp,
                                                  subspecies,
                                                  identifiedby,
                                                  latitude,
                                                  longitude)

diplusodon_coordClean <- diplusodon_coordClean %>% dplyr::select(id,
                                                institutioncode,
                                                collectioncode,
                                                catalognumber,
                                                gen_sp,
                                                subspecies,
                                                identifiedby,
                                                latitude,
                                                longitude)

diplusodon_clean <- rbind(diplusodon_coordClean, diplusodon_noCoord_inf)

rm(diplusodon_coordClean, diplusodon_noCoord_inf, i, j, correct.mun, generate.names, replace.names, titling)

#Writing *.csv 
#write.csv(diplusodon_clean, file = "Datasets/Diplusodon/diplusodon_clean.csv", row.names = FALSE)

#=====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Reading clean dataset
diplusodon_clean <- read.csv(file = "datasets/Diplusodon/diplusodon_clean.csv", na.strings = c("", NA))

#Standardizing gen_sp column
diplusodon_clean$gen_sp <- gsub(" ", "_", diplusodon_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the cr shapefile and, then, with the grids 
coords <- diplusodon_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- diplusodon_clean
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
write.csv(coords, "datasets/Diplusodon/diplusodon_cr.csv", row.names = FALSE)
