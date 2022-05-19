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

#Reading GBIF data (7,970)
trimezieae_gbif <- fread(file = "datasets/Trimezieae/0155773-200613084148143/occurrence.txt",
                        na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Reducing data dimensionality by selecting only necessary columns
trimezieae_gbif <- trimezieae_gbif %>% dplyr::select(institutionCode,
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
trimezieae_gbif <- trimezieae_gbif %>% rename("species" = specificEpithet,
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
trimezieae_gbif <- cbind(id = 1:nrow(trimezieae_gbif), trimezieae_gbif)

#=============#
# speciesLink #
#=============#

#Reading spLink (3,564)
trimezieae_spLink <- fread(file = "datasets/Trimezieae/speciesLink_all_110178_20210114214719.txt", 
                          na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Selecting important attributes
trimezieae_spLink <- trimezieae_spLink %>% dplyr::select(institutioncode,
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

#Coercing coords into numeric values (Warning:'NAs introduced by coercion' is expected because of non-numeric strings)
trimezieae_spLink$longitude <- as.numeric(as.character(trimezieae_spLink$longitude))
trimezieae_spLink$latitude <- as.numeric(as.character(trimezieae_spLink$latitude))

#Giving an unique ID number for each record
trimezieae_spLink <- cbind(id = (nrow(trimezieae_gbif) + 1):(nrow(trimezieae_gbif) + nrow(trimezieae_spLink)), trimezieae_spLink)

#=====================================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (27,280), and adding a column to define the original dataset
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

trimezieae <- merge.with.source(x = trimezieae_gbif,
                               y = trimezieae_spLink,
                               name.x = "gbif",
                               name.y = "splink")

rm(merge.with.source, trimezieae_gbif, trimezieae_spLink)

#=====================================================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (observations with NA for any of those attributes
#were not considered) (6,255)
a <- trimezieae
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
trimezieae <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
trimezieae <- trimezieae %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
trimezieae$latitude[trimezieae$latitude == 0] <- NA
trimezieae$longitude[trimezieae$longitude == 0] <- NA

#Removing registers without determiner's name (3,620) 
plyr::count(trimezieae$identifiedby)
trimezieae$identifiedby[trimezieae$identifiedby %in% c("?", "-", "#", "##", "//05/1987",
                                                     "24/X/2013", "`") | trimezieae$identifiedby == 0] <- NA
trimezieae <- trimezieae %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- trimezieae$stateprovince
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
trimezieae$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
trimezieae$stateprovince[trimezieae$stateprovince == "?" | trimezieae$stateprovince == "-"] <- NA
plyr::count(trimezieae$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur (2,596)
#according to Silveira et al (2016) 
trimezieae <- trimezieae %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                              "Minas Gerais", 
                                                                              "Bahia", 
                                                                              "Pernambuco", 
                                                                              "Paraiba", 
                                                                              "Mato Grosso",
                                                                              "Distrito Federal"))

#Removing records without species level identification (2,432)
trimezieae <- trimezieae %>% filter(!is.na(species))

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (2,409)
plyr::count(trimezieae$basisofrecord)
trimezieae <- trimezieae %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                                       "s", "S", "PreservedSpecimen", "PreserverdSpecimen"))
trimezieae <- trimezieae %>% dplyr::select(-basisofrecord)

#=====================================================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

#Creating a vector to work with
identifiedby <- trimezieae$identifiedby

#Extracting a vector including determiners' names. It is preferable to use this vector instead of the data set itself because any changes are easily reversible. 
identifiedby <- as.character(trimezieae$identifiedby)

#Ideally, it is better to start with a list of taxonomic specialists' compiled beforehand. Alternatively, as experts likely indentified the majority of samples from a given taxon, it is possible to infer specialists based on identification frequency. In this example we looked for specialists at the top of the list below. 
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#To improve accuracy, we confirmed if names in 'names.count' with at least 3 identifications were specialists by searching for taxonomic publications for the family of the focal group and authored by each name at Google Scholar. In this example, we searched for: allintitle: Iridaceae OR Trimezieae OR Neomarica OR Pseudotrimezia OR Trimezia OR Pseudiris author:"determiner".

#Next, based on the function 'replace.names', we standardized specialist's name. This is done in two iterations:
#(1) The first iteration returns, for manual evaluation, the automatically replaced names (names above the 'top' threshold) and names that are worth to be checked (names above the 'bottom' threshold but below the 'top' threshold).
#(2) In the second iteration, names that were erroneously replaced in the first iteration should be included in the argument 'not replace'. Likewise, names that were supposed to be replaced but were below the 'top' threshold should be included in the argument 'replace'.

#Because the procedure is iterative, the user should not change the vector 'identifiedBy' directly in the first iteration. For this purpose, we used a secondary vector ('identifiedBy_2'). After the second iteration, when everything should be set, the user must assign 'identifiedBy' as 'identifiedBy_2' before running the protocol for the following name. 

#At the end of the standardizing procedure, the following vector should contain all specialists' names
specialists <- c()

#NS Chukr
replace.by <- "NS Chukr"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("N. S. Churk","Nadio Saud Chukr",
                                            "Nadia Said Chukr",
                                            "Nadio Saud Chuckr",
                                            "Naudio Saud Chuckr"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#PF Ravenna
replace.by <- "PF Ravenna"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Ravenna/Publ: Bol.Soc.ARgBot, 10:321.1965",
                                            "R. Ravenna"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#J Lovo
replace.by <- "J Lovo"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("L. Lovo", "Lovo (2014)","J. Lou[...] 2015-03-12",
                                            "J. Lovo (SPF)","Mello-Silva, R; Lovo, J"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#ASB Gil
replace.by <- "ASB Gil"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("A. Gil/24-V-2007","A. Gil confirmação por J. Lovo.",
                                            "A. Gil/IX-2008", "A. S. B. Gil","André Gil (por foto)",
                                            "A. Gil/25-V-2007","A.Gil","A. dos S. Bragança Gil",
                                            "Ferreira Junior , C.A. confirmada por Gil, A."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#PN Oliveira
replace.by <- "PN Oliveira"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Oliveira, MS","Oliveira, RS","Oliveira, PA",
                                                "Oliveira, O.F. de"),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#L Capellari Júnior
replace.by <- "L Capellari Júnior"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("L. Capellari","L.Capellari, J.r.",
                                            "Lindolpho Capellari Jr","L.; Capellari Jr., L" ,
                                            "Lindolpho Capellari Jr."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#EBA Dias
replace.by <- "EBA Dias"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("E. B. A. Dias"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#MV Dantas-Queiroz
replace.by <- "MV Dantas-Queiroz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Queiroz, MD","; Queiroz, M",
                                            "M.Queiroz"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#L Eggers
replace.by <- "L Eggers"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#H Huaylla
replace.by <- "H Huaylla"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JE Henrich
replace.by <- "JE Henrich"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("J. E. Henrich; P. Goldblatt","J.E. Henrich & P. Goldblatt"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#B Mathew
replace.by <- "B Mathew"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#Replacing column
trimezieae$identifiedby <- identifiedby

#Filtering (1,945)
trimezieae <- trimezieae %>% filter(identifiedby %in% specialists)

rm(identifiedby, identifiedby_2, specialists, names.count, replace.by)

#======================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

library(flora)

#Generating a column for scientific names (without authors and including infraspecific epithet)
trimezieae$gen_sp <- paste(trimezieae$genus,
                          trimezieae$species,
                          trimezieae$subspecies,
                          sep = " ")

#Names
taxa <- plyr::count(trimezieae$gen_sp)
taxa <- as.character(taxa$x)

#Removing NA (character derived from 'subspecies' attribute)
for(i in 1:length(taxa)){
  taxa[i] <- gsub("NA", "", taxa[i])
  taxa[i] <- trimws(taxa[i])
}

#Suggesting with flora (and retrieving a few additional information that may be useful)
taxa_suggested <-get.taxa(taxa, vegetation.type = TRUE, 
                          habitat = TRUE, domain = TRUE, life.form = TRUE)

#Writing *.csv for manual checking
write.csv(taxa_suggested, file = "lists/Trimezieae/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Trimezieae/taxa_corrected.csv", stringsAsFactors = F, 
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

#Removing invalid taxa (1,821)
trimezieae$gen_sp <- gsub("NA", "", trimezieae$gen_sp)
trimezieae$gen_sp <- trimws(trimezieae$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen) | !taxa_gensp$gen %in% c("Neomarica",
                                                                                  "Trimezia",
                                                                                  "Pseudotrimezia",
                                                                                  "Pseudiris")]
trimezieae <- trimezieae[!trimezieae$gen_sp %in% invalid_taxa, ]

#Correcting the data set (varieties suppressed)
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(trimezieae)){
  for(j in 1:nrow(taxa_gensp)){
    if(trimezieae$gen_sp[i] == taxa_gensp$replace[j]){
      trimezieae$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

rm(taxa, taxa_suggested, taxa_corrected, taxa_gensp, invalid_taxa, str)

#======================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (849)
trimezieae_coord <- trimezieae %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (972)
trimezieae_noCoord <- trimezieae %>% filter(is.na(latitude) | is.na(longitude))

#Removing invalid coordinates and transferring them to trimezieae_noCoord (invalid coordinates are indicated by function clean_coordinates() in the next session)
trimezieae_noCoord <- rbind(trimezieae_noCoord, trimezieae_coord[c(701, 702, 703, 705, 706, 707), ])
trimezieae_coord <- trimezieae_coord[-c(701, 702, 703, 705, 706, 707), ]

#======================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (1,667)
trimezieae_coordFlagged <- trimezieae_coord %>% clean_coordinates(lon = "longitude",
                                                                lat = "latitude",
                                                                species = "gen_sp",
                                                                value = "flagged",
                                                                tests = c("equal", "gbif", 
                                                                          "institutions", 
                                                                          "outliers", "seas",
                                                                          "zeros"))

invalid_coords <- trimezieae_coord[trimezieae_coordFlagged == FALSE, ]
trimezieae_coordClean <- trimezieae_coord[trimezieae_coordFlagged  == TRUE, ]

#Binding invalid_coords to trimezieae_noCoord (4,511)
trimezieae_noCoord <- rbind(trimezieae_noCoord, invalid_coords)

rm(invalid_coords, trimezieae_coordFlagged, trimezieae_coord)

#======================================================================================#

#===============================#
# CLEANING MUNICIPALITIES NAMES #
#===============================#

#Dividing dataset in registers with and without values for stateprovince (72 and 937)
trimezieae_noState <- trimezieae_noCoord %>% filter(is.na(stateprovince))
trimezieae_state <- trimezieae_noCoord %>% filter(!is.na(stateprovince))

#NA values for municipalities (284)
trimezieae_na <- trimezieae_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(trimezieae_noState$county))
#names.count <- names.count[order(-names.count$freq), ]

#Reading list of municipalities according to DTB (Divisão Territorial Brasileira - IBGE) - ftp://geoftp.ibge.gov.br/organizacao_do_territorio/estrutura_territorial/divisao_territorial/
dtb <- fread(file = "lists/RELATORIO_DTB_BRASIL_MUNICIPIO.csv",
             na.strings = c("", NA), encoding = "UTF-8")
dtb <- dtb %>% filter(Nome_UF %in% c("Bahia", 
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
mun <- readOGR("shapefiles/br_municipios/BRMUE250GC_SIR.shp", encoding = "UTF-8", use_iconv = T) #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#Projecting
proj4string(cr) <- crswgs84
mun <- spTransform(mun, crswgs84)

#Defining municipalities in which the campos rupestres occur
list_mun <- over(cr, mun)
list_mun <- unique(list_mun$NM_MUNICIP)
list_mun <- as.character(list_mun) 

#Which municipalities are homonyms?
homonyms <- mun$NM_MUNICIP[mun$NM_MUNICIP %in% list_mun]
homonyms <- unique(correct.mun(homonyms[duplicated(homonyms)]))

#Standardizing municipalities names - lists and data sets
list_mun_std <- correct.mun(list_mun)
dtb$Nome_Município <- correct.mun(dtb$Nome_Município)
trimezieae_noState$municipality_gbif_std <- correct.mun(trimezieae_noState$municipality_gbif)
trimezieae_state$municipality_gbif_std <- correct.mun(trimezieae_state$municipality_gbif)
trimezieae_noState$county_std <- correct.mun(trimezieae_noState$county)
trimezieae_state$county_std <- correct.mun(trimezieae_state$county)

#Selecting dtb municipalities in which the campos rupestres occur
dtb_cr <- dtb[dtb$Nome_Município %in% list_mun_std, ]

#Manually checking and correcting by crossing the list assigned as 'check' with the dtb_cr dataset. 
check <- data.frame(plyr::count(trimezieae_noState$municipality_gbif_std)) 
trimezieae_noState$municipality_gbif_std[trimezieae_noState$municipality_gbif_std == 
                                          "Anaje"] <- "Anage"
trimezieae_noState$municipality_gbif_std[trimezieae_noState$municipality_gbif_std == 
                                          "Gouvea"] <- "Gouveia"
trimezieae_noState$municipality_gbif_std[trimezieae_noState$municipality_gbif_std == 
                                          "Sao Joao D'Alianca"] <- "Sao Joao D'alianca"
trimezieae_noState$municipality_gbif_std[trimezieae_noState$municipality_gbif_std 
                                        %in% c("Brasil", "Goias",
                                               "Mun?")] <- NA

check <- data.frame(plyr::count(trimezieae_noState$county_std)) 
trimezieae_noState$county_std[trimezieae_noState$county_std 
                             %in% c("Brasil", "5Km Sul da Cidade")] <- NA
trimezieae_na <- rbind(trimezieae_na, 
                      trimezieae_noState[which(is.na(trimezieae_noState$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                                is.na(trimezieae_noState$county_std)), 1:18]) 
trimezieae_noState <- trimezieae_noState[-which(is.na(trimezieae_noState$municipality_gbif_std) & #removing registers with NA values for these two attributes from trimezieae_noState
                                                is.na(trimezieae_noState$county_std)), ]

check <- data.frame(plyr::count(trimezieae_state$municipality_gbif_std))
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std %in% 
                                        c("Alto Garca")] <- "Alto Garcas"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std %in% 
                                        c("Alto Paraiso","Alto Paraiso Goias",
                                          "N de Alto Paraiso")] <- "Alto Paraiso de Goias"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Abaira   Piata"] <- "Abaira"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Anaje"] <- "Anage"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Barra de Mendes"] <- "Barra do Mendes"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Betania/Floresta"] <- "Betania"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std %in% 
                                        c("Braslilia")] <- "Brasilia"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Caitite"] <- "Caetite"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Catugi"] <- "Catuji"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Cajazeiras Sao Jose das Piranhas"] <- "Cajazeiras"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Cidade Ecletica"] <- "Santo Antonio do Descoberto"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std %in% 
                                        c("Conceicao do Mato de Dentro",
                                          "Conc do Mato Dentro")] <- "Conceicao do Mato Dentro"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std %in% 
                                        c("Conselheiro da Mata")] <- "Conselheiro Mata"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std %in% 
                                        c("Contenda do Sincora")] <- "Contendas do Sincora"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Corumba  de Goias"] <- "Corumba de Goias"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std %in% 
                                        c("Delfinopolis (?)","Delfinopoilis")] <- "Delfinopolis"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std %in% 
                                        c("Don Basilio")] <- "Dom Basilio"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Felixandia"] <- "Felixlandia"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Gouvea"] <- "Gouveia"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std %in% 
                                        c("Jaboticatuba")] <- "Jaboticatubas"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Olhos D#?#Agua"] <- "Olhos D'agua"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Planaltina de Goias"] <- "Planaltina"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std 
                                      %in% c("Presidente Kubitchek")] <- "Presidente Kubitschek"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Rio das Contas"] <- "Rio de Contas"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Santana do Pirapama"] <- "Santana de Pirapama"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std 
                                       %in% c("Santa Terezina")] <- "Santa Terezinha"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Sao Tome das Letras"] <- "Sao Thome das Letras"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Sao Joao da Alianca",
                                        "Sao Joao Dalianca"] <- "Sao Joao D'alianca"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std 
                                      %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                             "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std == 
                                        "Terezina de Goias", 
                                      "Mteresina de Goias"] <- "Teresina de Goias"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std 
                                      %in% c("Varzea de Palma")] <- "Varzea da Palma"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std 
                                      %in% c("Vila Bela Santissima Trindade", "Vila Bela da Sma Trindade",
                                             "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std 
                                      %in% c("Chapada da Cotagem", "Chapada Diamantina",
                                             "Chapada do Araripe", "Chapada dos Guimaraes",
                                             "Chapada dos Veadeiros", "Paraiba", "Pernambuco",
                                             "Serra do Cipo", "Serra do Espinhaco", "Serra do Tombador",
                                             "Ufba Ondina", "Veadeiros", "Bahia", "Chapadao do Ceu", "Goias", 
                                             "Minas Gerais", "Serra do Cipo", "Serra do Espinhaco", 
                                             "Serra do Tombador", "Veadeiros", "Bahia", "Betim/Brumadinho",
                                             "Coletada A Margem da Estrada de Macambinha","Br 135 Km 404",
                                             "Lavras Sao Joao Del Rey", "Ouro Preto/ Mariana",
                                             "Mato Grosso", "", "Chapada dos Veados", 
                                             "Ss Paraiso")] <- NA


check <- data.frame(plyr::count(trimezieae_state$county_std))
trimezieae_state$county_std[trimezieae_state$county_std %in% 
                             c("Alto Paraiso", "Alto Paraa-so de Goias",
                               "N de Alto Paraiso")] <- "Alto Paraiso de Goias"
trimezieae_state$county_std[trimezieae_state$county_std == 
                             "Anaje"] <- "Anage"
trimezieae_state$county_std[trimezieae_state$county_std == 
                             "Barra de Mendes"] <- "Barra do Mendes"
trimezieae_state$county_std[trimezieae_state$county_std == 
                             "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
trimezieae_state$county_std[trimezieae_state$county_std %in% 
                             c("Conceicao do Mato de Dentro",
                               "Conc do Mato Dentro","Conceicao Doato Dentro")] <- "Conceicao do Mato Dentro"
trimezieae_state$county_std[trimezieae_state$county_std == 
                             "Cristalina Mun"] <- "Cristalina"
trimezieae_state$county_std[trimezieae_state$county_std == 
                             "Don Basilio"] <- "Dom Basilio"
trimezieae_state$county_std[trimezieae_state$county_std %in% 
                             c("Golvea")] <- "Golveia"
trimezieae_state$county_std[trimezieae_state$county_std %in% 
                             c("Grao Mogol Mun")] <- "Grao Mogol"
trimezieae_state$county_std[trimezieae_state$county_std %in% 
                             c("Jaboticatuba")] <- "Jaboticatubas"
trimezieae_state$county_std[trimezieae_state$county_std == 
                             "Joaquim Fela-cio"] <- "Joaquim Felicio"
trimezieae_state$county_std[trimezieae_state$county_std %in% 
                             c("Morro do Chapeu Mun")] <- "Morro do Chapeu"
trimezieae_state$county_std[trimezieae_state$county_std %in% 
                             c("Planaltina de Goias")] <- "Planaltina"
trimezieae_state$county_std[trimezieae_state$county_std == 
                             "Rio das Contas",
                           "Rio de Contas Mun"] <- "Rio de Contas"
trimezieae_state$county_std[trimezieae_state$county_std %in% 
                             c("Sao Joao da Alianca", 
                               "Sao Joao D##Alianca")] <- "Sao Joao D'alianca"
trimezieae_state$county_std[trimezieae_state$county_std %in% 
                             c("Sao Tome das Letras ")] <- "Sao Thome das Letras "
trimezieae_state$county_std[trimezieae_state$county_std %in% 
                             c("Vila Bela de Santa-ssima Trindade", 
                               "Vila Bela de Santissima Trindade",
                               "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
trimezieae_state$county_std[trimezieae_state$county_std == 
                             "Santana de Pirapama"] <- "Santana do Pirapama"
trimezieae_state$municipality_gbif_std[trimezieae_state$municipality_gbif_std 
                                      %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                             "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
trimezieae_state$county_std[trimezieae_state$county_std == 
                             "Brasa-lia"] <- "Brasilia"
trimezieae_state$county_std[trimezieae_state$county_std 
                           %in% c("", "Chapada dos Guimaraes",
                                  "Chapada dos Veadeiros", "Chapada Gaucha",
                                  "Goias", "Minas Gerais", "No Disponible",
                                  "Pe", "Ufba Ondina", "Veadeiros", 
                                  "Minas Gerais", "Mato Grosso", "Serra de Ibitipoca",
                                  "Serra do Cabral", "Serra do Cipo", "Serra do Espinhaco",
                                  "Entre Serro E Lagoa Santa","Lavras Sao Joao Del Rey",
                                  "Chapada dos Veados")] <- NA

#Binding registers with NA values for these two attributes (2,791)
trimezieae_na <- rbind(trimezieae_na, 
                      trimezieae_state[which(is.na(trimezieae_state$municipality_gbif_std) & 
                                              is.na(trimezieae_state$county_std)), 1:18])

#Removing registers with NA values for these two attributes from trimezieae_state (2,438)
trimezieae_state <- trimezieae_state[which(!is.na(trimezieae_state$municipality_gbif_std) & 
                                           is.na(trimezieae_state$county_std)), ]

#Filtering trimezieae_noState and trimezieae_state with dtb_cr (29 and 1,762)
trimezieae_noState_filt <- trimezieae_noState %>% filter (municipality_gbif_std %in% dtb_cr$Nome_Município |
                                                          county_std %in% dtb_cr$Nome_Município)
trimezieae_state_filt <- trimezieae_state %>% filter (municipality_gbif_std %in% dtb_cr$Nome_Município |
                                                      county_std %in% dtb_cr$Nome_Município)

#=====================================================================================================#

#================================#
# INFERRING MUNICIPALITIES NAMES # 
#================================#

#Standardizing 'locality'
trimezieae_na$locality_std <- correct.mun(trimezieae_na$locality)

#Vector with municipalities names from dtb to be used with 'grepl'. This is important to match the whole municipality's name
grepl_munc_concla_cr <- c()
for(i in 1:nrow(dtb_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", dtb_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(trimezieae_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(trimezieae_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = trimezieae_na$locality_std[j])){
      vec[j] <- dtb_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(trimezieae_na$locality_std)] <- NA

#New column with inferred municipality names
trimezieae_na$municipality <- vec

#Removing registers with NA for 'municipality' (182)
trimezieae_na <- trimezieae_na %>% filter(!is.na(municipality))

#Concatenating information on municipality into a unique column 
trimezieae_noState_filt$municipality <- NA
for(i in 1:nrow(trimezieae_noState)){
  if(!is.na(trimezieae_noState_filt$municipality_gbif_std[i])){
    trimezieae_noState_filt$municipality[i] <- trimezieae_noState_filt$municipality_gbif_std[i] 
  } else if(!is.na(trimezieae_noState_filt$county_std[i])){
    trimezieae_noState_filt$municipality[i] <- trimezieae_noState_filt$county_std[i]
  }
}

trimezieae_state_filt$municipality <- NA
for(i in 1:nrow(trimezieae_state)){
  if(!is.na(trimezieae_state_filt$municipality_gbif_std[i])){
    trimezieae_state_filt$municipality[i] <- trimezieae_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(trimezieae_state_filt$county_std[i])){
    trimezieae_state_filt$municipality[i] <- trimezieae_state_filt$county_std[i]
  }
}

#Concatenating data sets without coordinates (484)
trimezieae_noState_filt <- trimezieae_noState_filt %>% 
  dplyr::select(colnames(trimezieae_noState_filt)[!colnames(trimezieae_noState_filt) %in%
                                                   c("county","county_std",
                                                     "municipality_gbif",
                                                     "municipality_gbif_std")])

trimezieae_state_filt <- trimezieae_state_filt %>% 
  dplyr::select(colnames(trimezieae_state_filt)[!colnames(trimezieae_state_filt) %in%
                                                 c("county","county_std",
                                                   "municipality_gbif",
                                                   "municipality_gbif_std")])

trimezieae_na <- trimezieae_na %>% 
  dplyr::select(colnames(trimezieae_na)[!colnames(trimezieae_na) %in%
                                         c("county","county_std",
                                           "municipality_gbif",
                                           "municipality_gbif_std")])

trimezieae_noCoord_inf <- rbind(trimezieae_noState_filt, trimezieae_state_filt, trimezieae_na, fill = TRUE)

#Registers occurring in homonyms municipalities (8)
reg_hom <- trimezieae_noCoord_inf[trimezieae_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
trimezieae_noCoord_inf <- trimezieae_noCoord_inf %>% filter(!municipality %in% homonyms)
trimezieae_noCoord_inf <- rbind(trimezieae_noCoord_inf, reg_hom) #binding them after removing those without state/province information

rm(trimezieae_na, trimezieae_noState, trimezieae_noState_filt, trimezieae_state, trimezieae_state_filt,
   dtb, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, trimezieae_noCoord, reg_hom, check)

#=====================================================================================================#

#===============================================#
# INFERRING COORDINATES BASED ON 'municipality' #
#===============================================#

#Loading centroids
mun_centr <- fread(file = "lists/Neotropics_municipalities.csv", na.strings = c("", NA))

#Filtering for Brazilian municipalities
mun_centr_br <- mun_centr %>% filter(country == "Brazil")

#Standardizing municipalities names
mun_centr_br$`minor area` <- correct.mun(mun_centr_br$`minor area`)

#Filtering for campos rupestres municipalities
mun_centr_cr <- mun_centr_br %>% filter(`minor area` %in% dtb_cr$Nome_Município &
                                          `major area` %in% dtb_cr$Nome_UF)

#Removing homonyms from the data set and replacing them 
centr_hom <- mun_centr_br %>% filter(`minor area` %in% homonyms)
mun_centr_cr <- mun_centr_cr %>% filter(!`minor area` %in% homonyms)

#Some municipalities (according to the DTB) are not included in centroids data set
not_included <- dtb_cr %>% filter(!Nome_Município %in% mun_centr_cr$`minor area`)

#Inferring coordinates based on 'municipality'
for(i in 1:length(trimezieae_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(trimezieae_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      trimezieae_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      trimezieae_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Inferring coordinates for homonyms municipalities. Here, I relied on the combination between state and municipality. If an specific record lacks information on any of those attributes, coordinates are not inferred.
for(i in 1:length(trimezieae_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(trimezieae_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       trimezieae_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(trimezieae_noCoord_inf$stateprovince[i])){
      trimezieae_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      trimezieae_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing records without coordinates (476)
trimezieae_noCoord_inf <- trimezieae_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, dtb_cr, centr_hom, homonyms, not_included)

#=================================================================================================#

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unimportant columns (1,288)
trimezieae_noCoord_inf <- trimezieae_noCoord_inf %>% dplyr::select(id,
                                                                 institutioncode,
                                                                 collectioncode,
                                                                 catalognumber,
                                                                 gen_sp,
                                                                 subspecies,
                                                                 identifiedby,
                                                                 latitude,
                                                                 longitude)

trimezieae_coordClean <- trimezieae_coordClean %>% dplyr::select(id,
                                                               institutioncode,
                                                               collectioncode,
                                                               catalognumber,
                                                               gen_sp,
                                                               subspecies,
                                                               identifiedby,
                                                               latitude,
                                                               longitude)

trimezieae_clean <- rbind(trimezieae_coordClean, trimezieae_noCoord_inf)

#Is there any duplicated registers? If so, removing them (1,206)
trimezieae_clean[duplicated(trimezieae_clean$id) == TRUE, ]
trimezieae_clean <- unique(trimezieae_clean) 

rm(trimezieae, trimezieae_coordClean, trimezieae_noCoord_inf, i, j, correct.mun, generate.names,
   replace.names, titling)

#====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Standardizing gen_sp column
trimezieae_clean$gen_sp <- gsub(" ", "_", trimezieae_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the campos rupestres shapefile and, then, with the grids (1,686)
coords <- trimezieae_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- trimezieae_clean
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
coords <- coords[ , -c("id_2")]

#Saving dataset
write.csv(coords, "datasets/Trimezieae/trimezieae_cr.csv", row.names = FALSE)
