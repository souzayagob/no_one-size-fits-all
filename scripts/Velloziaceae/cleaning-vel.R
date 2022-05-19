library(tidyverse) 
library(data.table) 
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 

source('C:/Users/Yago/Dropbox/Mestrado/Bancos_organizing/Scripts/functions.R') #Windows
source('~/Dropbox/Mestrado/Bancos_organizing/Scripts/functions.R') #Ubuntu
#======================================================================================#

#==================#
# READING DATASETS #
#==================#

#======#
# GBIF #
#======#

#Reading gbif data
vel_gbif <- fread(file = "Datasets/Velloziaceae/0029865-200613084148143-gbif_velloziaceae/occurrence.txt",
                     na.strings = c("", NA), stringsAsFactors = FALSE, encoding = "UTF-8")

#Selecting important attributes
vel_gbif <- vel_gbif %>% select(institutionCode,
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
vel_gbif <- vel_gbif %>% rename("species" = specificEpithet,
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
vel_gbif <- cbind(id = 1:nrow(vel_gbif), vel_gbif)

#Encoding strings as UTF-8 (R version 4)
#class(vel_gbif) <- "data.frame" #if the current class is "data.table" "data.frame",
#the following code does not work properly. 
#for(i in 1:ncol(vel_gbif)){
  #if(is.character(vel_gbif[ , i])){
 #   Encoding(vel_gbif[ , i]) <- "UTF-8" 
#  }
#}

#=============#
# speciesLink #
#=============#

#Reading spLink 
vel_spLink <- fread(file = "Datasets/Velloziaceae/speciesLink_all_105625_20200728173230.txt", 
                       na.strings = c("", NA), stringsAsFactors = FALSE, encoding = "UTF-8")

#Selecting important attributes
vel_spLink <- vel_spLink %>% select(institutioncode,
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
vel_spLink$longitude <- as.numeric(as.character(vel_spLink$longitude))
vel_spLink$latitude <- as.numeric(as.character(vel_spLink$latitude))

#Giving an unique ID number for each record
vel_spLink <- cbind(id = (nrow(vel_gbif) + 1):(nrow(vel_gbif) + nrow(vel_spLink)), vel_spLink)

#Encoding strings as UTF-8 
#class(vel_spLink) <- "data.frame" #if the current class is "data.table" "data.frame",
#the following code does not work properly. 
#for(i in 1:ncol(vel_spLink)){
  #if(is.character(vel_spLink[ , i])){
   # Encoding(vel_spLink[ , i]) <- "UTF-8" 
  #}
#}

#==============================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (28,519), and adding a column to define the original dataset
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

vel <- merge.with.source(x = vel_gbif,
                            y = vel_spLink,
                            name.x = "gbif",
                            name.y = "splink")

rm(merge.with.source, vel_gbif, vel_spLink)

#==============================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (observations with NA for any of those attributes
#were not considered) (22,562). Apparently, there is no duplicated register. Is that right? 
a <- vel
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
vel <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
vel <- vel %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
vel$latitude[vel$latitude == 0] <- NA
vel$longitude[vel$longitude == 0] <- NA

#Removing registers without identifier name (10,922)
#plyr::count(vel$identifiedby)
vel$identifiedby[vel$identifiedby %in% c("?", "-", "#", "##", "//05/1987",
                                               "24/X/2013") | vel$identifiedby == 0] <- NA
vel <- vel %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- vel$stateprovince
#lookup_states <- data.frame("Incorreto" = plyr::count(province)[ , 1])
lookup_states <- fread(file = "lists/lookup_states.csv", na.strings = c("", NA))
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
vel$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
vel$stateprovince[vel$stateprovince == "?" | vel$stateprovince == "-"] <- NA
plyr::count(vel$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur 
#according to Silveira et al (2016) (9,849). NA's are not considered.
vel <- vel %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                        "Minas Gerais", 
                                                                        "Bahia", 
                                                                        "Pernambuco", 
                                                                        "Paraiba", 
                                                                        "Mato Grosso",
                                                                        "Distrito Federal"))

#Removing records without species level identification (9,133)
vel <- vel %>% filter(!is.na(species))
vel <- vel %>% filter(!species %in% c("sp.", "sp1"))
#plyr::count(vel$species)

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (9,078)
#plyr::count(vel$basisofrecord)
vel <- vel %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                                 "s", "S", "PreservedSpecimen"))
vel <- vel %>% select(-basisofrecord)

#==============================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

#Creating a vector to work with
identifiedby <- vel$identifiedby

#Generating a list of identifiers ordered by frequency of identifications (more identifications
#is supposed to be related with a better identification)
#Counting check
#names.count <- as.data.frame(plyr::count(identifiedby))
#names.count <- names.count[order(-names.count$freq), ]

#R Mello-Silva
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                              check.by = generate.names("R.", "Mello-Silva"),
                              replace.by = "R Mello-Silva",
                              not.replace = c("Mello, TRB"),
                              replace = c("R, Mello-Silva confirmada","R, Mello -Silva,confirmada por",
                                          "R, Mello-Silva Confirmada por","R, Mello-Silva confirmada por",
                                          "R, Mello-Silva,confirmada por","R, Mello-Silva,Confirmada por",
                                          "R, Mello -Silva, Confirmada por","Renato Mello -Silva",
                                          "R. Melo Silva; N. L. Menezes","Renato de Mello-Silva",
                                          "R, Mello-Silva, confirmada por","R, Mello -Silva confirmada",
                                          "Renato Mello-Silva","Silva, R. M.",
                                          "R.M. Silva & M. Menezes","Renato Mello-Silva Confirmada por",
                                          "Renato Mello-Silva confirmada por",
                                          "R.M. Silva","Menezes & Mello-Silva","Silva, R.M.","Zappi, DC; Mello-Silva, R",
                                          "Renato Mello -Silva confirmada por","Renato de Mello-Silva; b.1961; Mello-Silva",
                                          "fide Mello-Sila & Sasaki","Silva, R.M. da","R.M. Silva & N. Menezes",
                                          "Silva, RM","fide Mello-Silva & Sasaki","F.Mello-Silva (SPF)",
                                          "Silva, RM da"))

#N Menezes
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("N.", "Menezes"),
                              replace.by = "N Menezes",
                              not.replace = c(),
                              replace = c("N. L. de Menezes","Silva, R.M.; Menezes, N.L.",
                                          "N. L. de Menezes; J. Semir","M. L. Menezes",
                                          "N. L. Menezes; R. Mello Silva"))

#LB Smith
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                              check.by = generate.names("L.B.", "Smith"),
                              replace.by = "LB Smith",
                              not.replace = c("Smith, S.F.","Smith, LD"),
                              replace = c("Smith, LB; Ayensu, ES","L. B. Smith & E.S Ayensu",
                                          "Smith & Ayensu","L.B.Smith et Ayensu",
                                          "Smith; Ayensu, LB","L. B. Smith; E. S, Ayensu",
                                          "Ayensu; Smith, L.B.","L. B. Smith & E. S. Ayensu",
                                          "LBS","\"\"Ayensu, E.S.; Smith, L.B.\"\"",
                                          "L. B. Smith; Ayensu","L. B. Smith; E. S. Ayensu",
                                          "L. B. Smith & Ayensu","L. B. Smith; E. Ayensu",
                                          "Ayensu, E; Smith, LB","Ayensu, ES; Smith, LB"))

#RJV Alves
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                              check.by = generate.names("R.J.V.", "Alves"),
                              replace.by = "RJV Alves",
                              not.replace = c("Alves, M.V."),
                              replace = c())

#L Montserrat
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                              check.by = generate.names("L.", "Montserrat"),
                              replace.by = "L Montserrat",
                              not.replace = c(),
                              replace = c())

#Goethart
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                              check.by = c("Goethart", "GOETHART"),
                              replace.by = "Goethart",
                              not.replace = c(),
                              replace = c())

#write.csv(names.count, file = "identificadores_vel.csv", row.names = FALSE)

#Replacing column
vel$identifiedby <- identifiedby

#Filtering (7,044)
specialists <- c("R Mello-Silva",
                 "N Menezes",
                 "LB Smith",
                 "RJV Alves",
                 "L Montserrat",
                 "Goethart")
vel <- vel %>% filter(identifiedby %in% specialists)

rm(identifiedby, specialists)

#====================================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

#Generating a column for scientific names (without authors and including infraspecific epithet)
vel$gen_sp <- paste(vel$genus,
                       vel$species,
                       vel$subspecies,
                       sep = " ")

#taxa <- plyr::count(vel$gen_sp)
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
#write.csv(taxa_suggested, file = "lists/Velloziaceae/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Velloziaceae/taxa_corrected.csv", stringsAsFactors = F, 
                           na.strings = c("NA",""))


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

#Removing invalid taxa (7,006)
vel$gen_sp <- gsub("NA", "", vel$gen_sp)
vel$gen_sp <- trimws(vel$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen)]
vel <- vel[!vel$gen_sp %in% invalid_taxa, ]

#Correcting the dataset
taxa_gensp$replace <- taxa_corrected$taxa
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(vel)){
  for(j in 1:nrow(taxa_gensp)){
    if(vel$gen_sp[i] == taxa_gensp$replace[j]){
      vel$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

#Checking if there is any last corrections to be made
vel$gen_sp[!vel$gen_sp %in% paste(taxa_gensp$gen, taxa_gensp$sp, sep = " ")]

rm(taxa_corrected, taxa_gensp, invalid_taxa, str)

#======================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (2,731)
vel_coord <- vel %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (4,275)
vel_noCoord <- vel %>% filter(is.na(latitude) | is.na(longitude))

#======================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (2625)
vel_coordFlagged <- vel_coord %>% clean_coordinates(lon = "longitude",
                                                          lat = "latitude",
                                                          species = "gen_sp",
                                                          value = "flagged",
                                                          tests = c("equal", "gbif", 
                                                                    "institutions", 
                                                                    "outliers", "seas",
                                                                    "zeros"))

invalid_coords <- vel_coord[vel_coordFlagged == FALSE, ]
vel_coordClean <- vel_coord[vel_coordFlagged  == TRUE, ]

#Binding invalid_coords to vel_noCoord (4,381)
vel_noCoord <- rbind(vel_noCoord, invalid_coords)

rm(invalid_coords, vel_coordFlagged, vel_coord)

#Writing *.csv 
#write.csv(vel_coordClean, file = "Datasets/Velloziaceae/vel_coordClean.csv", row.names = FALSE)
#======================================================================================#

#===============================#
# CLEANING MUNICIPALITIES NAMES #
#===============================#

#Dividing dataset in registers with and without values for stateprovince (494 and 3,887)
vel_noState <- vel_noCoord %>% filter(is.na(stateprovince))
vel_state <- vel_noCoord %>% filter(!is.na(stateprovince))

#NA values for municipalities (1,670)
vel_na <- vel_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(vel_noState$county))
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
vel_noState$municipality_gbif_std <- correct.mun(vel_noState$municipality_gbif)
vel_state$municipality_gbif_std <- correct.mun(vel_state$municipality_gbif)
vel_noState$county_std <- correct.mun(vel_noState$county)
vel_state$county_std <- correct.mun(vel_state$county)

#Selecting mun_concla municipalities in which the campos rupestres occur
mun_concla_cr <- mun_concla[mun_concla$Nome_Município %in% list_mun_std, ]

#nrow(mun_concla_cr) is not equal to length(list_mun_std). Why?
#Municipalities from the São Paulo state
#list_mun_std[!list_mun_std %in% mun_concla_cr$Nome_Município]

#Manually checking and correcting by crossing the list assigned as 'check' with the mun_concla_cr dataset. 
#check <- data.frame(plyr::count(vel_noState$municipality_gbif_std)) #ok
vel_noState$municipality_gbif_std[vel_noState$municipality_gbif_std == 
                                    "Anaje"] <- "Anage"
vel_noState$municipality_gbif_std[vel_noState$municipality_gbif_std == 
                                       "Gouvea"] <- "Gouveia"
vel_noState$municipality_gbif_std[vel_noState$municipality_gbif_std == 
                                     "Sao Joao D'Alianca"] <- "Sao Joao D'alianca"
vel_noState$municipality_gbif_std[vel_noState$municipality_gbif_std 
                                %in% c("Brasil", "Goias")] <- NA

#check <- data.frame(plyr::count(vel_noState$county_std)) #ok
vel_noState$county_std[vel_noState$county_std 
                                  %in% c("Brasil", "5Km Sul da Cidade")] <- NA
vel_na <- rbind(vel_na, 
                vel_noState[which(is.na(vel_noState$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                    is.na(vel_noState$county_std)), 1:18]) 
vel_noState <- vel_noState[-which(is.na(vel_noState$municipality_gbif_std) & #removing registers with NA values for these two attributes from vel_noState
                                  is.na(vel_noState$county_std)), ]

#check <- data.frame(plyr::count(vel_state$municipality_gbif_std))
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std %in% 
                                     c("Alto Paraiso","Alto Paraiso Goias")] <- "Alto Paraiso de Goias"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                  "Abaira   Piata"] <- "Abaira"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                     "Anaje"] <- "Anage"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                     "Barra de Mendes"] <- "Barra do Mendes"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                     "Betania/Floresta"] <- "Betania"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                     "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                  "Catugi"] <- "Catuji"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                  "Cajazeiras Sao Jose das Piranhas"] <- "Cajazeiras"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                  "Cidade Ecletica"] <- "Santo Antonio do Descoberto"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                  "Conceicao do Mato de Dentro"] <- "Conceicao do Mato Dentro"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                     "Corumba  de Goias"] <- "Corumba de Goias"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                  "Felixandia"] <- "Felixlandia"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                     "Gouvea"] <- "Gouveia"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                  "Olhos D#?#Agua"] <- "Olhos D'agua"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                     "Planaltina de Goias"] <- "Planaltina"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                     "Rio das Contas"] <- "Rio de Contas"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                  "Santana do Pirapama"] <- "Santana de Pirapama"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                     "Sao Tome das Letras"] <- "Sao Thome das Letras"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                     "Sao Joao da Alianca"] <- "Sao Joao D'alianca"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std 
                                %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei")] <- "Sao Joao Del Rei"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std == 
                                  "Terezina de Goias"] <- "Teresina de Goias"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std 
                                %in% c("Varzea de Palma")] <- "Varzea da Palma"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std 
                                %in% c("Vila Bela Santissima Trindade", "Vila Bela da Sma Trindade",
                                       "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
vel_state$municipality_gbif_std[vel_state$municipality_gbif_std 
                       %in% c("Chapada da Cotagem", "Chapada Diamantina",
                              "Chapada do Araripe", "Chapada dos Guimaraes",
                              "Chapada dos Veadeiros", "Paraiba", "Pernambuco",
                              "Serra do Cipo", "Serra do Espinhaco", "Serra do Tombador",
                              "Ufba Ondina", "Veadeiros", "Bahia", "Chapadao do Ceu", "Goias", 
                              "Minas Gerais", "Serra do Cipo", "Serra do Espinhaco", 
                              "Serra do Tombador", "Veadeiros", "Bahia", "Betim/Brumadinho",
                              "Coletada A Margem da Estrada de Macambinha")] <- NA


#check <- data.frame(plyr::count(vel_state$county_std))
vel_state$county_std[vel_state$county_std %in% 
                          c("Alto Paraiso", "Alto Paraa-so de Goias")] <- "Alto Paraiso de Goias"
vel_state$county_std[vel_state$county_std == 
                          "Anaje"] <- "Anage"
vel_state$county_std[vel_state$county_std == 
                          "Barra de Mendes"] <- "Barra do Mendes"
vel_state$county_std[vel_state$county_std == 
                          "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
vel_state$county_std[vel_state$county_std == 
                          "Cristalina Mun"] <- "Cristalina"
vel_state$county_std[vel_state$county_std %in% 
                       c("Grao Mogol Mun")] <- "Grao Mogol"
vel_state$county_std[vel_state$county_std == 
                          "Joaquim Fela-cio"] <- "Joaquim Felicio"
vel_state$county_std[vel_state$county_std %in% 
                       c("Morro do Chapeu Mun")] <- "Morro do Chapeu"
vel_state$county_std[vel_state$county_std == 
                          "Rio das Contas"] <- "Rio de Contas"
vel_state$county_std[vel_state$county_std %in% 
                          c("Sao Joao da Alianca", 
                            "Sao Joao D##Alianca")] <- "Sao Joao D'alianca"
vel_state$county_std[vel_state$county_std %in% 
                       c("Sao Tome das Letras ")] <- "Sao Thome das Letras "
vel_state$county_std[vel_state$county_std %in% 
                          c("Vila Bela de Santa-ssima Trindade", 
                            "Vila Bela de Santissima Trindade",
                            "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
vel_state$county_std[vel_state$county_std == 
                          "Santana de Pirapama"] <- "Santana do Pirapama"
vel_state$county_std[vel_state$county_std == 
                          "Brasa-lia"] <- "Brasilia"
vel_state$county_std[vel_state$county_std 
                                %in% c("", "Chapada dos Guimaraes",
                                       "Chapada dos Veadeiros", "Chapada Gaucha",
                                       "Goias", "Minas Gerais", "No Disponible",
                                       "Pe", "Ufba Ondina", "Veadeiros", 
                                       "Minas Gerais", "Mato Grosso", "Serra de Ibitipoca",
                                       "Serra do Cabral", "Serra do Cipo", "Serra do Espinhaco")] <- NA
vel_na <- rbind(vel_na, 
                vel_state[which(is.na(vel_state$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                    is.na(vel_state$county_std)), 1:18]) 
vel_state <- vel_state[-which(is.na(vel_state$municipality_gbif_std) & #removing registers with NA values for these two attributes from vel_noState
                                    is.na(vel_state$county_std)), ]

#Filtering vel_noState and vel_state with mun_concla_cr
vel_noState_filt <- vel_noState %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                                    county_std %in% mun_concla_cr$Nome_Município)
vel_state_filt <- vel_state %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                                county_std %in% mun_concla_cr$Nome_Município)

#====================================================================================================#

#================================#
# INFERRING MUNICIPALITIES NAMES #
#================================#

#Standardizing 'locality'
vel_na$locality_std <- correct.mun(vel_na$locality)

#Vector with municipalities names from mun_concla to be used with grepl
grepl_munc_concla_cr <- c()
for(i in 1:nrow(mun_concla_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", mun_concla_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(vel_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(vel_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = vel_na$locality_std[j])){
      vec[j] <- mun_concla_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(vel_na$locality_std)] <- NA

#New column with inferred municipality names
vel_na$municipality <- vec

#Removing registers with NA for 'municipality'
vel_na <- vel_na %>% filter(!is.na(municipality))

#Concatenating information on municipality into a unique column 
#vel_noState_filt$municipality <- NA
#for(i in 1:nrow(vel_noState)){
 # if(!is.na(vel_noState_filt$municipality_gbif_std[i])){
  #  vel_noState_filt$municipality[i] <- vel_noState_filt$municipality_gbif_std[i] 
  #} else if(!is.na(vel_noState_filt$county_std[i])){
  #  vel_noState_filt$municipality[i] <- vel_noState_filt$county_std[i]
  #}
#}

vel_state_filt$municipality <- NA
for(i in 1:nrow(vel_state)){
  if(!is.na(vel_state_filt$municipality_gbif_std[i])){
    vel_state_filt$municipality[i] <- vel_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(vel_state_filt$county_std[i])){
    vel_state_filt$municipality[i] <- vel_state_filt$county_std[i]
  }
}

#Concatenating data sets
#vel_noState_filt <- vel_noState_filt[ , -which(names(vel_noState_filt) 
                                                     #%in% c("county","county_std",
                                                      #      "municipality_gbif",
                                                       #     "municipality_gbif_std"))]

vel_state_filt <- vel_state_filt[ , -which(names(vel_state_filt) 
                                                 %in% c("county","county_std",
                                                        "municipality_gbif",
                                                        "municipality_gbif_std"))]
vel_na <- vel_na[ , -which(names(vel_na) 
                                 %in% c("county", "municipality_gbif", "locality_std"))]

vel_noCoord_inf<- rbind(vel_noState_filt, vel_state_filt, vel_na)

#Registers occurring in homonyms municipalities
reg_hom <- vel_noCoord_inf[vel_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
vel_noCoord_inf <- vel_noCoord_inf %>% filter(!municipality %in% homonyms)
vel_noCoord_inf <- rbind(vel_noCoord_inf, reg_hom) #binding them after removing those without state/province information

rm(vel_na, vel_noState, vel_noState_filt, vel_state, vel_state_filt,
   mun_concla, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, vel_noCoord, reg_hom)

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
for(i in 1:length(vel_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(vel_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      vel_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      vel_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Homonyms
for(i in 1:length(vel_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(vel_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       vel_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(vel_noCoord_inf$stateprovince[i])){
      vel_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      vel_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing NA's for coordinates
vel_noCoord_inf <- vel_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, mun_concla_cr, centr_hom, homonyms)

#=================================================================================================#

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unnimportant columns
vel_noCoord_inf <- vel_noCoord_inf %>% select(id,
                                                    institutioncode,
                                                    collectioncode,
                                                    catalognumber,
                                                    gen_sp,
                                                    subspecies,
                                                    identifiedby,
                                                    latitude,
                                                    longitude)

vel_coordClean <- vel_coordClean %>% select(id,
                                                  institutioncode,
                                                  collectioncode,
                                                  catalognumber,
                                                  gen_sp,
                                                  subspecies,
                                                  identifiedby,
                                                  latitude,
                                                  longitude)

vel_clean <- rbind(vel_coordClean, vel_noCoord_inf)

#Is there any duplicated registers? 
#vel_clean[duplicated(vel_clean$id) == TRUE, ] #no

rm(vel, vel_coordClean, vel_noCoord_inf, i, j, correct.mun, generate.names,
   replace.names, titling)

#Writing *.csv 
#write.csv(vel_clean, file = "Datasets/Velloziaceae/vel_clean.csv", row.names = FALSE)

#=====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Reading clean dataset
vel_clean <- read.csv(file = "datasets/Velloziaceae/vel_clean.csv", na.strings = c("", NA))

#Standardizing gen_sp column
vel_clean$gen_sp <- gsub(" ", "_", vel_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the cr shapefile and, then, with the grids 
coords <- vel_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- vel_clean
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
write.csv(coords, "datasets/Velloziaceae/vel_cr.csv", row.names = FALSE)
