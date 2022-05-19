#====================================================================================================#

library(tidyverse) 
library(data.table) 
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 
library(raster)

#Loading functions
source('scripts/functions.R') 

#====================================================================================================#

#==================#
# READING DATASETS #
#==================#

#======#
# GBIF #
#======#

#Reading GBIF data (2,920)
minaria_gbif <- fread(file = "datasets/Minaria/0155691-200613084148143/occurrence.txt",
                     na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Reducing data dimensionality by selecting only necessary columns
minaria_gbif <- minaria_gbif %>% dplyr::select(institutionCode,
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

#Renaming attributes in order to match speciesLink column names
minaria_gbif <- minaria_gbif %>% rename("species" = specificEpithet,
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
minaria_gbif <- cbind(id = 1:nrow(minaria_gbif), minaria_gbif)

#=============#
# speciesLink #
#=============#

#Reading spLink (1,661)
minaria_spLink <- fread(file = "datasets/Minaria/speciesLink_all_65777_20210114200802.txt", 
                       na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Selecting important attributes
minaria_spLink <- minaria_spLink %>% dplyr::select(institutioncode,
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

#Coercing coords into numeric values 
minaria_spLink$longitude <- as.numeric(as.character(minaria_spLink$longitude))
minaria_spLink$latitude <- as.numeric(as.character(minaria_spLink$latitude))

#Giving an unique ID number for each record
minaria_spLink <- cbind(id = (nrow(minaria_gbif) + 1):(nrow(minaria_gbif) + nrow(minaria_spLink)), minaria_spLink)

#======================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (4,581), and adding a column to define the original dataset
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

minaria <- merge.with.source(x = minaria_gbif,
                            y = minaria_spLink,
                            name.x = "gbif",
                            name.y = "splink")

rm(merge.with.source, minaria_gbif, minaria_spLink)

#======================================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (observations with NA for any of those attributes
#were not considered) (3,614)
a <- minaria
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
minaria <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
minaria <- minaria %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
minaria$latitude[minaria$latitude == 0] <- NA
minaria$longitude[minaria$longitude == 0] <- NA

#Removing registers without identifier name (2,531) 
#plyr::count(minaria$identifiedby)
minaria$identifiedby[minaria$identifiedby %in% c("?", "-", "#", "##", "//05/1987",
                                               "24/X/2013") | minaria$identifiedby == 0] <- NA
minaria <- minaria %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- minaria$stateprovince
#lookup_states <- data.frame("Incorreto" = plyr::count(province)[ , 1])
lookup_states <- fread(file = "lists/lookup_states.csv", na.strings = c("", NA),
                       encoding = "UTF-8")
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
minaria$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
minaria$stateprovince[minaria$stateprovince == "?" | minaria$stateprovince == "-"] <- NA
plyr::count(minaria$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur (2,273)
#according to Silveira et al (2016) 
minaria <- minaria %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                        "Minas Gerais", 
                                                                        "Bahia", 
                                                                        "Pernambuco", 
                                                                        "Paraiba", 
                                                                        "Mato Grosso",
                                                                        "Distrito Federal"))

#Removing records without species level identification (2,270)
minaria <- minaria %>% filter(!is.na(species))

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (2,269)
#plyr::count(minaria$basisofrecord)
minaria <- minaria %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                                 "s", "S", "PreservedSpecimen"))
minaria <- minaria %>% dplyr::select(-basisofrecord)

#=====================================================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

#Creating a vector to work with
identifiedby <- minaria$identifiedby

#Extracting a vector including determiners' names. It is preferable to use this vector instead of the data set itself because any changes are easily reversible. 
identifiedby <- as.character(minaria$identifiedby)

#Ideally, it is better to start with a list of taxonomic specialists' compiled beforehand. Alternatively, as experts likely indentified the majority of samples from a given taxon, it is possible to infer specialists based on identification frequency. In this example we looked for specialists at the top of the list below. 
names.count <- as.data.frame(plyr::count(identifiedby))
names.count <- names.count[order(-names.count$freq), ]

#To improve accuracy, we confirmed if names in 'names.count' with at least 3 identifications were specialists by searching for taxonomic publications for the family of the focal group and authored by each name at Google Scholar. In this example, we searched for: Apocynaceae OR Minaria author:"determiner".

#Next, based on the function 'replace.names', we standardized specialist's name. This is done in two iterations:
#(1) The first iteration returns, for manual evaluation, the automatically replaced names (names above the 'top' threshold) and names that are worth to be checked (names above the 'bottom' threshold but below the 'top' threshold).
#(2) In the second iteration, names that were erroneously replaced in the first iteration should be included in the argument 'not replace'. Likewise, names that were supposed to be replaced but were below the 'top' threshold should be included in the argument 'replace'.

#Because the procedure is iterative, the user should not change the vector 'identifiedBy' directly in the first iteration. For this purpose, we used a secondary vector ('identifiedBy_2'). After the second iteration, when everything should be set, the user must assign 'identifiedBy' as 'identifiedBy_2' before running the protocol for the following name. 

#A Rapini
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("A", "Rapini"),
                                replace.by = "A Rapini",
                                not.replace = c(),
                                replace = c("Alessandro Rapini (por foto)","A. Rapini 2000-03",
                                            "A. Rapini/10-06-2010", "Fernandes , M.G.C, confirmada por A. Rapini",
                                            "Fernandes , M.G.C, confirmada por A. Rapini",
                                            "Fernandes, M.G.C. confirmada por Rapini,A.","Fernandes, MG.C confirmada por A. Rapini",
                                            "Fernandes, MGC, confirmada por Rapini,A."))
identifiedby <- identifiedby_2

#S Liede-Schumann
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("S", "Liede-Schumann"),
                                replace.by = "S Liede-Schumann",
                                not.replace = c(),
                                replace = c("S. Liede; U. Meve", "S. Liede"))
identifiedby <- identifiedby_2

#WD Stevens
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("WD", "Stevens"),
                                replace.by = "WD Stevens",
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2

#MA Farinaccio
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("MA", "Farinaccio"),
                                replace.by = "MA Farinaccio",
                                not.replace = c(),
                                replace = c("Maria Ana Farinaccio, VII."))
identifiedby <- identifiedby_2

#J Fontella-Pereira
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("J", "Fontella-Pereira"),
                                replace.by = "J Fontella-Pereira",
                                not.replace = c(),
                                replace = c("Fontella J.","J. Fontella","J.FONTELLA",
                                            "G. Fontella","J. Fontella, 9/II/","J. Fontella 1983",
                                            "J. Fontella 1991","J. Fontella 2005-07-09",
                                            "J.P.Fontella & N.F.da S.Marquete/VI-1975",
                                            "J.Fontella","P.J. Fontella","J. Fontella P.",
                                            "Fontellah 1979-08-20","J.P. Fontella 1978-05-25",
                                            "J. Fontella 2005-07","Fontella, J.P.; Marquete, N."))
identifiedby <- identifiedby_2

#TUP Konno
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("TUP", "Konno"),
                                replace.by = "TUP Konno",
                                not.replace = c(),
                                replace = c("TATIANA KONNO","Tatiana Konno","Tatiana konno"))
identifiedby <- identifiedby_2

#C Bitencourt
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("C", "Bitencourt"),
                                replace.by = "C Bitencourt",
                                not.replace = c(),
                                replace = c("C.BITENCOURT"))
identifiedby <- identifiedby_2

#PL Ribeiro
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("PL", "Ribeiro"),
                                replace.by = "PL Ribeiro",
                                not.replace = c(),
                                replace = c("Patrícia Luz Ribeiro"))
identifiedby <- identifiedby_2

#DJ Goyder
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("DJ", "Goyder"),
                                replace.by = "DJ Goyder",
                                not.replace = c(),
                                replace = c("Dr David Goyder", "D. Goyder 1997",
                                            "G. D. Goyder 1997"))
identifiedby <- identifiedby_2

#AO Simões
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("AO", "Simões" ),
                                replace.by = "AO Simões",
                                not.replace = c(),
                                replace = c("Simão, A O"))
identifiedby <- identifiedby_2

#JF Pereira
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("JF", "Pereira"),
                                replace.by = "JF Pereira",
                                not.replace = c("Pereira, L.C."),
                                replace = c())
identifiedby <- identifiedby_2

#RGP Santos
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names("RGP", "Santos"),
                                replace.by = "RGP Santos",
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2

#Replacing column
minaria$identifiedby <- identifiedby

#Filtering (2,046)
specialists <- c("A Rapini",
                 "S Liede-Schumann",
                 "WD Stevens",
                 "MA Farinaccio",
                 "J Fontella-Pereira",
                 "TUP Konno",
                 "C Bitencourt",
                 "PL Ribeiro",
                 "DJ Goyder",
                 "AO Simões",
                 "JF Pereira",
                 "RGP Santos")
minaria <- minaria %>% filter(identifiedby %in% specialists)

rm(identifiedby, identifiedby_2, specialists, names.count)

#======================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

library(flora)

#Generating a column for scientific names (without authors and including infraspecific epithet)
minaria$gen_sp <- paste(minaria$genus,
                       minaria$species,
                       minaria$subspecies,
                       sep = " ")

#Names
taxa <- plyr::count(minaria$gen_sp)
taxa <- as.character(taxa$x)

#Removing NA (character derived from 'subspecies' attribute)
for(i in 1:length(taxa)){
 taxa[i] <- gsub("NA", "", taxa[i])
 taxa[i] <- trimws(taxa[i])
}

#Suggesting with flora (and retrieving a few additional informations that may be useful)
taxa_suggested <-get.taxa(taxa, vegetation.type = TRUE, 
                          habitat = TRUE, domain = TRUE, life.form = TRUE)

#Writing *.csv for manual checking
#write.csv(taxa_suggested, file = "lists/Minaria/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Minaria/taxa_corrected.csv", stringsAsFactors = F, 
                           na.strings = c("NA",""))

#Establishing a data set with information on genus, species and varieties
taxa_gensp <- tibble(gen = NA, sp = NA, var = NA, .rows = nrow(taxa_corrected))
for(i in 1:nrow(taxa_corrected)){
  str <- strsplit(taxa_corrected$corrected.name[i], split = " ")[[1]]
  if(length(str) == 2){
    taxa_gensp$gen[i] <- str[1]
    taxa_gensp$sp[i] <- str[2]
  } else if(length(str) == 4){
    taxa_gensp$gen[i] <- str[1]
    taxa_gensp$sp[i] <- str[2]
    taxa_gensp$var[i] <- paste(str[3], str[4], sep = " ")
  }
}

#Original names and correspondent corrected names 
taxa_gensp$replace <- taxa_corrected$original.search

#Removing invalid taxa (2,039)
minaria$gen_sp <- gsub("NA", "", minaria$gen_sp)
minaria$gen_sp <- trimws(minaria$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen) | taxa_gensp$gen != "Minaria"]
minaria <- minaria[!minaria$gen_sp %in% invalid_taxa, ]

#Correcting the data set
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(minaria)){
  for(j in 1:nrow(taxa_gensp)){
    if(minaria$gen_sp[i] == taxa_gensp$replace[j]){
      minaria$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

rm(taxa, taxa_suggested, taxa_corrected, taxa_gensp, invalid_taxa, str)

#======================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (717)
minaria_coord <- minaria %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (1,322)
minaria_noCoord <- minaria %>% filter(is.na(latitude) | is.na(longitude))

#Removing invalid coordinates and transfering them to minaria_noCoord
minaria_noCoord <- rbind(minaria_noCoord, minaria_coord[c(562, 563, 564, 565,
                                                          566, 567, 568, 569), ])
minaria_coord <- minaria_coord[-c(562, 563, 564, 565,
                                  566, 567, 568, 569), ]


#======================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (684)
minaria_coordFlagged <- minaria_coord %>% clean_coordinates(lon = "longitude",
                                                          lat = "latitude",
                                                          species = "gen_sp",
                                                          value = "flagged",
                                                          tests = c("equal", "gbif", 
                                                                    "institutions", 
                                                                    "outliers", "seas",
                                                                    "zeros"))

invalid_coords <- minaria_coord[minaria_coordFlagged == FALSE, ]
minaria_coordClean <- minaria_coord[minaria_coordFlagged  == TRUE, ]

#Binding invalid_coords to minaria_noCoord (1,355)
minaria_noCoord <- rbind(minaria_noCoord, invalid_coords)

rm(invalid_coords, minaria_coordFlagged, minaria_coord)

#======================================================================================#

#===============================#
# CLEANING MUNICIPALITIES NAMES #
#===============================#

#Dividing dataset in registers with and without values for stateprovince (98 and 1,257)
minaria_noState <- minaria_noCoord %>% filter(is.na(stateprovince))
minaria_state <- minaria_noCoord %>% filter(!is.na(stateprovince))

#NA values for municipalities (343)
minaria_na <- minaria_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(minaria_noState$county))
#names.count <- names.count[order(-names.count$freq), ]

#Reading list of municipalities according to CTB (Divisão Territorial Brasileira - IBGE) - ftp://geoftp.ibge.gov.br/organizacao_do_territorio/estrutura_territorial/divisao_territorial/
mun_concla <- fread(file = "lists/RELATORIO_DTB_BRASIL_MUNICIPIO.csv",
                    na.strings = c("", NA))
mun_concla <- mun_concla %>% filter(Nome_UF %in% c("Bahia", 
                                                   "Goiás",
                                                   "Mato Grosso",
                                                   "Minas Gerais",
                                                   "Pernambuco",
                                                   "Paraíba",
                                                   "Distrito Federal"))

#Defining projection
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

#Reading shapefiles
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
mun <- readOGR("shapefiles/br_municipios/BRMUE250GC_SIR.shp") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#Projecting
proj4string(cr) <- crswgs84
mun <- spTransform(mun, crswgs84)

#Defining municipalities in which the campos rupestres occur
list_mun <- over(cr, mun)
list_mun <- unique(list_mun$NM_MUNICIP)
list_mun <- as.character(list_mun) 

#Which municipalities are homonyms?
homonyms <- mun$NM_MUNICIP[mun$NM_MUNICIP %in% list_mun]
homonyms <- correct.mun(homonyms[duplicated(homonyms)])

#Standardizing municipalities names - lists and data sets
list_mun_std <- correct.mun(list_mun)
mun_concla$Nome_Município <- correct.mun(mun_concla$Nome_Município)
minaria_noState$municipality_gbif_std <- correct.mun(minaria_noState$municipality_gbif)
minaria_state$municipality_gbif_std <- correct.mun(minaria_state$municipality_gbif)
minaria_noState$county_std <- correct.mun(minaria_noState$county)
minaria_state$county_std <- correct.mun(minaria_state$county)

#Selecting mun_concla municipalities in which the campos rupestres occur
mun_concla_cr <- mun_concla[mun_concla$Nome_Município %in% list_mun_std, ]

#Manually checking and correcting by crossing the list assigned as 'check' with the mun_concla_cr dataset. 
#check <- data.frame(plyr::count(minaria_noState$municipality_gbif_std)) #ok
minaria_noState$municipality_gbif_std[minaria_noState$municipality_gbif_std == 
                                     "Anaje"] <- "Anage"
minaria_noState$municipality_gbif_std[minaria_noState$municipality_gbif_std == 
                                     "Gouvea"] <- "Gouveia"
minaria_noState$municipality_gbif_std[minaria_noState$municipality_gbif_std == 
                                     "Sao Joao D'Alianca"] <- "Sao Joao D'alianca"
minaria_noState$municipality_gbif_std[minaria_noState$municipality_gbif_std 
                                   %in% c("Brasil", "Goias",
                                          "Mun?")] <- NA

#check <- data.frame(plyr::count(minaria_noState$county_std)) #ok
minaria_noState$county_std[minaria_noState$county_std 
                        %in% c("Brasil", "5Km Sul da Cidade")] <- NA
minaria_na <- rbind(minaria_na, 
                 minaria_noState[which(is.na(minaria_noState$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                      is.na(minaria_noState$county_std)), 1:18]) 
minaria_noState <- minaria_noState[-which(is.na(minaria_noState$municipality_gbif_std) & #removing registers with NA values for these two attributes from minaria_noState
                                      is.na(minaria_noState$county_std)), ]

#check <- data.frame(plyr::count(minaria_state$municipality_gbif_std))
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std %in% 
                                   c("Alto Garca")] <- "Alto Garcas"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std %in% 
                                   c("Alto Paraiso","Alto Paraiso Goias")] <- "Alto Paraiso de Goias"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Abaira   Piata"] <- "Abaira"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Anaje"] <- "Anage"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Barra de Mendes"] <- "Barra do Mendes"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Betania/Floresta"] <- "Betania"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std %in% 
                                   c("Braslilia")] <- "Brasilia"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Catugi"] <- "Catuji"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Cajazeiras Sao Jose das Piranhas"] <- "Cajazeiras"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Cidade Ecletica"] <- "Santo Antonio do Descoberto"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std %in% 
                                   c("Conceicao do Mato de Dentro",
                                     "Conc do Mato Dentro")] <- "Conceicao do Mato Dentro"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Corumba  de Goias"] <- "Corumba de Goias"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std %in% 
                                   c("Delfinopolis (?)","Delfinopoilis")] <- "Delfinopolis"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Felixandia"] <- "Felixlandia"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Gouvea"] <- "Gouveia"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std %in% 
                                   c("Jaboticatuba")] <- "Jaboticatubas"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Olhos D#?#Agua"] <- "Olhos D'agua"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Planaltina de Goias"] <- "Planaltina"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std 
                                 %in% c("Presidente Kubitchek")] <- "Presidente Kubitschek"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Rio das Contas"] <- "Rio de Contas"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Santana do Pirapama"] <- "Santana de Pirapama"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Sao Tome das Letras"] <- "Sao Thome das Letras"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Sao Joao da Alianca"] <- "Sao Joao D'alianca"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std 
                                 %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                        "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std == 
                                   "Terezina de Goias"] <- "Teresina de Goias"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std 
                                 %in% c("Varzea de Palma")] <- "Varzea da Palma"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std 
                                 %in% c("Vila Bela Santissima Trindade", "Vila Bela da Sma Trindade",
                                        "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std 
                                 %in% c("Chapada da Cotagem", "Chapada Diamantina",
                                        "Chapada do Araripe", "Chapada dos Guimaraes",
                                        "Chapada dos Veadeiros", "Paraiba", "Pernambuco",
                                        "Serra do Cipo", "Serra do Espinhaco", "Serra do Tombador",
                                        "Ufba Ondina", "Veadeiros", "Bahia", "Chapadao do Ceu", "Goias", 
                                        "Minas Gerais", "Serra do Cipo", "Serra do Espinhaco", 
                                        "Serra do Tombador", "Veadeiros", "Bahia", "Betim/Brumadinho",
                                        "Coletada A Margem da Estrada de Macambinha","Br 135 Km 404",
                                        "Lavras Sao Joao Del Rey", "Ouro Preto/ Mariana")] <- NA


#check <- data.frame(plyr::count(minaria_state$county_std))
minaria_state$county_std[minaria_state$county_std %in% 
                        c("Alto Paraiso", "Alto Paraa-so de Goias")] <- "Alto Paraiso de Goias"
minaria_state$county_std[minaria_state$county_std == 
                        "Anaje"] <- "Anage"
minaria_state$county_std[minaria_state$county_std == 
                        "Barra de Mendes"] <- "Barra do Mendes"
minaria_state$county_std[minaria_state$county_std == 
                        "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std %in% 
                                   c("Conceicao do Mato de Dentro",
                                     "Conc do Mato Dentro","Conceicao Doato Dentro")] <- "Conceicao do Mato Dentro"
minaria_state$county_std[minaria_state$county_std == 
                        "Cristalina Mun"] <- "Cristalina"
minaria_state$county_std[minaria_state$county_std %in% 
                        c("Golvea")] <- "Golveia"
minaria_state$county_std[minaria_state$county_std %in% 
                        c("Grao Mogol Mun")] <- "Grao Mogol"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std %in% 
                                   c("Jaboticatuba")] <- "Jaboticatubas"
minaria_state$county_std[minaria_state$county_std == 
                        "Joaquim Fela-cio"] <- "Joaquim Felicio"
minaria_state$county_std[minaria_state$county_std %in% 
                        c("Morro do Chapeu Mun")] <- "Morro do Chapeu"
minaria_state$county_std[minaria_state$county_std == 
                        "Rio das Contas"] <- "Rio de Contas"
minaria_state$county_std[minaria_state$county_std %in% 
                        c("Sao Joao da Alianca", 
                          "Sao Joao D##Alianca")] <- "Sao Joao D'alianca"
minaria_state$county_std[minaria_state$county_std %in% 
                        c("Sao Tome das Letras ")] <- "Sao Thome das Letras "
minaria_state$county_std[minaria_state$county_std %in% 
                        c("Vila Bela de Santa-ssima Trindade", 
                          "Vila Bela de Santissima Trindade",
                          "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
minaria_state$county_std[minaria_state$county_std == 
                        "Santana de Pirapama"] <- "Santana do Pirapama"
minaria_state$municipality_gbif_std[minaria_state$municipality_gbif_std 
                                 %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                        "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
minaria_state$county_std[minaria_state$county_std == 
                        "Brasa-lia"] <- "Brasilia"
minaria_state$county_std[minaria_state$county_std 
                      %in% c("", "Chapada dos Guimaraes",
                             "Chapada dos Veadeiros", "Chapada Gaucha",
                             "Goias", "Minas Gerais", "No Disponible",
                             "Pe", "Ufba Ondina", "Veadeiros", 
                             "Minas Gerais", "Mato Grosso", "Serra de Ibitipoca",
                             "Serra do Cabral", "Serra do Cipo", "Serra do Espinhaco",
                             "Entre Serro E Lagoa Santa","Lavras Sao Joao Del Rey")] <- NA

minaria_na <- rbind(minaria_na, 
                 minaria_state[which(is.na(minaria_state$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                    is.na(minaria_state$county_std)), 1:18]) 
minaria_state <- minaria_state[-which(is.na(minaria_state$municipality_gbif_std) & #removing registers with NA values for these two attributes from minaria_noState
                                  is.na(minaria_state$county_std)), ]

#Filtering minaria_noState and minaria_state with mun_concla_cr
minaria_noState_filt <- minaria_noState %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                                county_std %in% mun_concla_cr$Nome_Município)
minaria_state_filt <- minaria_state %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                            county_std %in% mun_concla_cr$Nome_Município)

#=====================================================================================================#

#================================#
# INFERRING MUNICIPALITIES NAMES # 
#================================#

#Standardizing 'locality'
minaria_na$locality_std <- correct.mun(minaria_na$locality)

#Vector with municipalities names from mun_concla to be used with grepl
grepl_munc_concla_cr <- c()
for(i in 1:nrow(mun_concla_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", mun_concla_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(minaria_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(minaria_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = minaria_na$locality_std[j])){
      vec[j] <- mun_concla_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(minaria_na$locality_std)] <- NA

#New column with inferred municipality names
minaria_na$municipality <- vec

#Removing registers with NA for 'municipality'
minaria_na <- minaria_na %>% filter(!is.na(municipality))

#Concatenating information on municipality into a unique column 
minaria_noState_filt$municipality <- NA
for(i in 1:nrow(minaria_noState)){
  if(!is.na(minaria_noState_filt$municipality_gbif_std[i])){
    minaria_noState_filt$municipality[i] <- minaria_noState_filt$municipality_gbif_std[i] 
  } else if(!is.na(minaria_noState_filt$county_std[i])){
    minaria_noState_filt$municipality[i] <- minaria_noState_filt$county_std[i]
  }
}

minaria_state_filt$municipality <- NA
for(i in 1:nrow(minaria_state)){
  if(!is.na(minaria_state_filt$municipality_gbif_std[i])){
    minaria_state_filt$municipality[i] <- minaria_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(minaria_state_filt$county_std[i])){
    minaria_state_filt$municipality[i] <- minaria_state_filt$county_std[i]
  }
}

#Concatenating data sets 
minaria_noState_filt <- minaria_noState_filt %>% 
  dplyr::select(colnames(minaria_noState_filt)[!colnames(minaria_noState_filt) %in%
                                                      c("county","county_std",
                                                        "municipality_gbif",
                                                        "municipality_gbif_std")])

minaria_state_filt <- minaria_state_filt %>% 
  dplyr::select(colnames(minaria_state_filt)[!colnames(minaria_state_filt) %in%
                                                 c("county","county_std",
                                                   "municipality_gbif",
                                                   "municipality_gbif_std")])

minaria_na <- minaria_na %>% 
  dplyr::select(colnames(minaria_na)[!colnames(minaria_na) %in%
                                                 c("county","county_std",
                                                   "municipality_gbif",
                                                   "municipality_gbif_std")])

minaria_noCoord_inf<- rbind(minaria_noState_filt, minaria_state_filt, minaria_na, fill = TRUE)

#Registers occurring in homonyms municipalities
reg_hom <- minaria_noCoord_inf[minaria_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
minaria_noCoord_inf <- minaria_noCoord_inf %>% filter(!municipality %in% homonyms)
minaria_noCoord_inf <- rbind(minaria_noCoord_inf, reg_hom) #binding them after removing those without state/province information

rm(minaria_na, minaria_noState, minaria_noState_filt, minaria_state, minaria_state_filt,
   mun_concla, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, minaria_noCoord, reg_hom)

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
for(i in 1:length(minaria_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(minaria_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      minaria_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      minaria_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Homonyms
for(i in 1:length(minaria_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(minaria_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       minaria_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(minaria_noCoord_inf$stateprovince[i])){
      minaria_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      minaria_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing NA's for coordinates (1,002)
minaria_noCoord_inf <- minaria_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, mun_concla_cr, centr_hom, homonyms)

#=================================================================================================#

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unnimportant columns (1,686)
minaria_noCoord_inf <- minaria_noCoord_inf %>% dplyr::select(id,
                                                institutioncode,
                                                collectioncode,
                                                catalognumber,
                                                gen_sp,
                                                subspecies,
                                                identifiedby,
                                                latitude,
                                                longitude)

minaria_coordClean <- minaria_coordClean %>% dplyr::select(id,
                                              institutioncode,
                                              collectioncode,
                                              catalognumber,
                                              gen_sp,
                                              subspecies,
                                              identifiedby,
                                              latitude,
                                              longitude)

minaria_clean <- rbind(minaria_coordClean, minaria_noCoord_inf)

#Is there any duplicated registers? If so, removing them (1,562)
minaria_clean[duplicated(minaria_clean$id) == TRUE, ]
minaria_clean <- unique(minaria_clean) 

rm(minaria, minaria_coordClean, minaria_noCoord_inf, i, j, correct.mun, generate.names,
   replace.names, titling)

#Writing *.csv 
write.csv(minaria_clean, file = "datasets/Minaria/minaria_clean.csv", row.names = FALSE)
#====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Reading clean dataset
minaria_clean <- read.csv(file = "datasets/Minaria/minaria_clean.csv", na.strings = c("", NA))

#Standardizing gen_sp column
minaria_clean$gen_sp <- gsub(" ", "_", minaria_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the cr shapefile and, then, with the grids 
coords <- minaria_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- minaria_clean
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
write.csv(coords, "datasets/Minaria/minaria_cr.csv", row.names = FALSE)
