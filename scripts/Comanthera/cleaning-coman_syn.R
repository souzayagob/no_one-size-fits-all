library(tidyverse) 
library(data.table) 
library(CoordinateCleaner) 
library(flora) 
library(rgdal) 

source('~/Dropbox/Mestrado/Bancos_organizing/Scripts/functions.R')

#======================================================================================#

#==================#
# READING DATASETS #
#==================#

#======#
# GBIF #
#======#

#Reading gbif data
coman_syn_gbif <- fread(file = "Datasets/Comanthera/0009147-191105090559680_coman_syn_gbif/occurrence.txt",
                    na.strings = c("", NA))

#Selecting important attributes
coman_syn_gbif <- coman_syn_gbif %>% dplyr::select(institutionCode,
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
coman_syn_gbif <- coman_syn_gbif %>% rename("species" = specificEpithet,
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
coman_syn_gbif <- cbind(id = 1:nrow(coman_syn_gbif), coman_syn_gbif)

#=============#
# speciesLink #
#=============#

#Reading spLink data
coman_syn_spLink <- fread(file = "Datasets/Comanthera/coman_syn_spLink.txt", 
                      na.strings = c("", NA))

#Selecting important attributes
coman_syn_spLink <- coman_syn_spLink %>% dplyr::select(institutioncode,
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
coman_syn_spLink$longitude <- as.numeric(as.character(coman_syn_spLink$longitude))
coman_syn_spLink$latitude <- as.numeric(as.character(coman_syn_spLink$latitude))

#Giving an unique ID number for each record
coman_syn_spLink <- cbind(id = (nrow(coman_syn_gbif) + 1):(nrow(coman_syn_gbif) + nrow(coman_syn_spLink)), coman_syn_spLink)

#======================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (24,278), and adding a column to define the original dataset
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

coman_syn <- merge.with.source(x = coman_syn_gbif,
                           y = coman_syn_spLink,
                           name.x = "gbif",
                           name.y = "splink")

rm(merge.with.source, coman_syn_gbif, coman_syn_spLink)

#======================================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (NA's not considered) (19,341)
a <- coman_syn
a.prime <- a[!is.na(a$institutioncode) &
               !is.na(a$collectioncode) &
               !is.na(a$catalognumber), ]
a.na <- a[is.na(a$institutioncode) |
            is.na(a$collectioncode) |
            is.na(a$catalognumber), ]
a <- a.prime[order(a.prime$genus), ]
a <- a[!duplicated(a, by = c("institutioncode", 
                            "collectioncode",
                            "catalognumber")), ]

a <- rbind(a, a.na)
coman_syn <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
coman_syn <- coman_syn %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
coman_syn$latitude[coman_syn$latitude == 0] <- NA
coman_syn$longitude[coman_syn$longitude == 0] <- NA

#Removing registers without identifier name (11,944)
#plyr::count(coman_syn$identifiedby)
coman_syn$identifiedby[coman_syn$identifiedby %in% c("?", "-") | coman_syn$identifiedby == 0] <- NA
coman_syn <- coman_syn %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- coman_syn$stateprovince
#lookup_states <- data.frame("Incorreto" = plyr::count(province)[ , 1], 
#                            "Correto" = NA)
lookup_states <- fread(file = "lists/lookup_states.csv", na.strings = c("", NA))
#write.csv(lookup_states, file = "lists/lookup_states_coman_syn.csv", row.names = F)
get_states <- lookup_states$Incorreto
names(get_states) <- lookup_states$Correto
for(i in 1:length(province)){
  for(j in 1:length(get_states)){
    if(province[i] == unname(get_states[j]) & !is.na(province[i])){
      province[i] <- names(get_states[j])
    }
  }
}
province <- trimws(province)
coman_syn$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
coman_syn$stateprovince[coman_syn$stateprovince == "?" | coman_syn$stateprovince == "-"] <- NA
#plyr::count(coman_syn$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur 
#according to Silveira et al (2016) (7,487)
coman_syn <- coman_syn %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                      "Minas Gerais", 
                                                                      "Bahia", 
                                                                      "Pernambuco", 
                                                                      "Paraiba", 
                                                                      "Mato Grosso",
                                                                      "Distrito Federal"))

#Removing records without species level identification (7,224)
#plyr::count(coman_syn$species)
coman_syn <- coman_syn %>% filter(!is.na(species))
coman_syn <- coman_syn %>% filter(!species %in% c("sp.", "sp1"))


#Removing the attribute 'basisofrecord' (7,222)
#plyr::count(coman_syn$basisofrecord) 
coman_syn <- coman_syn %>% filter(!is.na(basisofrecord))
coman_syn <- coman_syn %>% dplyr::select(-basisofrecord)

#======================================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

identifiedby <- coman_syn$identifiedby

#Generating a list of identifiers ordered by frequency of identifications (more identifications
#is supposed to be related with a better identification)
#Counting check
#names.count <- as.data.frame(plyr::count(identifiedby))
#names.count <- names.count[order(-names.count$freq), ]

#L. Echternacht
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("L.", "Echternacht"),
                              replace.by = "L. Echternacht",
                              not.replace = c(),
                              replace = c("Lívia Echternacht", "Livia Echternacht",
                                          "Livia A. Echternacht", "L.E. Andrande",
                                          "Andrade, E.L.", "L.Echternacht & M.Trovó",
                                          "Trovó, M.L.O.; Echternacht, L.", "M. Trovó; L. Echternacht"),
                              return = TRUE)

#M. T. C. Watanabe
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M. T. C.", "Watanabe"),
                              replace.by = "M. T. C. Watanabe",
                              not.replace = c(),
                              replace = c("M.Watanam", "M. Watanary",
                                          "M. Watamary", "M.Watanary",
                                          "M.Watamary", "N. Watamary",
                                          "M.WATANABE & M.TROVÓ", "M. Watanabe (SPF)",
                                          "Maurício Takashi Coutinho Watanabe",
                                          "Maurício Takahi Coutinho Watanabe"),
                              return = TRUE)

#A. M. Giulietti
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("A. M.", "Giulietti"),
                              replace.by = "A. M. Giulietti",
                              not.replace = c(),
                              replace = c("A.M. Giulietti & N. Hensold",
                                          "A.M.Giulietti Harley",
                                          "A.M.Giulietti & A.C.Pereira",
                                          "A.C.Pereira & A.M.Giulietti",
                                          "Ana Maria Giulietti",
                                          "Ana Maria Giulietti; b.1945; Giul.",
                                          "A.C.S.Pereira & A.M.Giulietti",
                                          "A.M. Giulietti & L.R. Parra",
                                          "A.M. Giulieet, L.R. Parra"),
                              return = T)

#L. R. Parra
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("L. R.", "Parra"),
                              replace.by = "L. R. Parra",
                              not.replace = c(),
                              replace = c(),
                              return = T)

#H. N. Moldenke
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("H. N", "Moldenke"),
                              replace.by = "H. N. Moldenke",
                              not.replace = c(),
                              replace = c("Harold N. Moldenke"),
                              return = T)

#B. Mourão
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("B.", "Mourão"),
                              replace.by = "B. Mourão",
                              not.replace = c(),
                              replace = c(),
                              return = T)

#M. Trovó
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M.", "Trovó"),
                              replace.by = "M. Trovó",
                              not.replace = c(),
                              replace = c("M.Trovó & C. Sarquis",
                                          "Trovó, MLO; Sano, PT; Echternacht, L",
                                          "M.L. Trovó & C. Sarquis",
                                          "M.L.O. Trovó/04-XI-2014",
                                          "M.L.O. Trovó; C. Sarquis",
                                          "C.Sarquis & M.Trovó",
                                          "Marcelo Trovó"),
                              return = T)

#P. T. Sano
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("P. T.", "Sano"),
                              replace.by = "P. T. Sano",
                              not.replace = c(),
                              replace = c("Sano, PT; Trovó, MLO",
                                          "Sano, P.T.; Costa, F.N.; Trovó, M.L.O.; Echternacht, L.; Andrino, C."),
                              return = T)

#N. Hensold
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("N.", "Hensold"),
                              replace.by = "N. Hensold",
                              not.replace = c(),
                              replace = c("Nancy Hensold", "hensold"),
                              return = T)

#L. J. Sauthier
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("L. J.", "Sauthier"),
                              replace.by = "L. J. Sauthier",
                              not.replace = c(),
                              replace = c(),
                              return = T)

#V. L. Scatena
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("V. L.", "Scatena"),
                              replace.by = "V. L. Scatena",
                              not.replace = c(),
                              replace = c(),
                              return = T)

#D. C. Zappi
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("D. C.", "Zappi"),
                              replace.by = "D. C. Zappi",
                              not.replace = c(),
                              replace = c(),
                              return = T)

#R. M. Harley
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("R. M.", "Harley"),
                              replace.by = "R. M. Harley",
                              not.replace = c(),
                              replace = c(),
                              return = T)

#R. C. Forzza
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("R. C.", "Forzza"),
                              replace.by = "R. C. Forzza",
                              not.replace = c(),
                              replace = c(),
                              return = T)

#A. Oliveira
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("A.", "Oliveira"),
                              replace.by = "A. Oliveira",
                              not.replace = c("Oliveira, G.B.",
                                              "Oliveira, M.", 
                                              "Oliveira, M",
                                              "Oliveira, MS", "Oliveira, RC",
                                              "Oliveira, MS; Chacon, RG",
                                              "Oliveira, G.C.", "Oriani, A.",
                                              "Oliveira, M.T.L.",
                                              "A. Silveira"),
                              replace = c(),
                              return = T)

#S. Splett 
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("S.", "Splett"),
                              replace.by = "S. Splett",
                              not.replace = c(),
                              replace = c(),
                              return = T)



#Counting check
#names.count <- as.data.frame(plyr::count(identifiedby))
#names.count <- names.count[order(-names.count$freq), ]

#Replacing column
coman_syn$identifiedby <- identifiedby

#Filtering (6,161)
specialists <- c("L. Echternacht",
                 "M. T. C. Watanabe",
                 "A. M. Giulietti",
                 "L. R. Parra",
                 "H. N. Moldenke",
                 "B. Mourão",
                 "M. Trovó",
                 "P. T. Sano",
                 "N. Hensold",
                 "L. J. Sauthier",
                 "V. L. Scatena",
                 "D. C. Zappi",
                 "R.M. Harley",
                 "R. C. Forzza",
                 "A. Oliveira",
                 "S. Splett")
coman_syn <- coman_syn %>% filter(identifiedby %in% specialists)

rm(identifiedby, specialists)

#====================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

#Generating a column for scientific names (without authors and including infraspecific epithet)
coman_syn$gen_sp <- paste(coman_syn$genus,
                      coman_syn$species,
                      coman_syn$subspecies,
                      sep = " ")

#taxa <- plyr::count(coman_syn$gen_sp)
#taxa$x <- as.character(taxa$x)
#taxa_suggested <- tibble("taxa" = taxa$x, "suggested" = NA)

#Removing NA (character derived from 'subspecies' attribute)
#for(i in 1:nrow(taxa_suggested)){
 # taxa_suggested[i, 1] <- gsub("NA", "", taxa_suggested[i, 1])
  #taxa_suggested[i, 1] <- trimws(taxa_suggested[i, 1])
#}

#Suggesting with flora (some errors)
#for(i in 1:nrow(taxa_suggested)){
 # taxa_suggested$suggested[i] <- suggest.names(taxa_suggested$taxa[i])
#}

#Writing *.csv for manual checking
#write.csv(taxa_suggested, file = "lists/Comanthera/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Comanthera/taxa_corrected.csv", stringsAsFactors = F, 
                           na.strings = c("NA",""))

#trimws for every column
for(i in 1:ncol(taxa_corrected)){
  taxa_corrected[ , i] <- trimws(taxa_corrected[ , i])
}

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

#Removing invalid taxa (6,161)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen)]
coman_syn <- coman_syn[!coman_syn$gen_sp %in% invalid_taxa, ]

#Correcting the dataset
coman_syn$gen_sp <- trimws(gsub("NA", "", coman_syn$gen_sp))
taxa_gensp$replace <- taxa_corrected$taxa
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(coman_syn)){
  for(j in 1:nrow(taxa_gensp)){
    if(coman_syn$gen_sp[i] == taxa_gensp$replace[j]){
      coman_syn$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

#Filtering for Comanthera (2,771)
coman <- coman_syn[which(grepl("Comanthera", coman_syn$gen_sp, fixed = TRUE)), ]

rm(taxa_corrected, taxa_gensp, str, coman_syn, invalid_taxa)
#=====================================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (1,375)
coman_coord <- coman %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (1,396)
coman_noCoord <- coman %>% filter(is.na(latitude) | is.na(longitude))

#=====================================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (1,267)
coman_coordFlagged <- coman_coord %>% clean_coordinates(lon = "longitude",
                                                          lat = "latitude",
                                                          species = "gen_sp",
                                                          value = "flagged",
                                                          tests = c("equal", "gbif", 
                                                                    "institutions", 
                                                                    "outliers", "seas",
                                                                    "zeros"))


invalid_coords <- coman_coord[coman_coordFlagged == FALSE, ]
coman_coordClean <- coman_coord[coman_coordFlagged  == TRUE, ]

#Binding invalid_coords to coman_noCoord
coman_noCoord <- rbind(coman_noCoord, invalid_coords)

rm(invalid_coords, coman_coordFlagged, coman_coord)

#Writing *.csv
write.csv(coman_coordClean, file = "Datasets/Comanthera/coman_coordClean.csv", row.names = FALSE)
#====================================================================================================#

#===============================#
# CLEANING MUNICIPALITIES NAMES #
#===============================#

#Dividing dataset in registers with and without values for stateprovince (280 and 1,224)
coman_noState <- coman_noCoord %>% filter(is.na(stateprovince))
coman_state <- coman_noCoord %>% filter(!is.na(stateprovince))

#NA values for municipalities (447)
coman_na <- coman_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(coman_noState$county))
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

#Testing intersection
#plot(mun[mun$NM_MUNICIP %in% list_mun, ]) 

#Which municipalities are homonyms?
homonyms <- mun$NM_MUNICIP[mun$NM_MUNICIP %in% list_mun]
homonyms <- correct.mun(homonyms[duplicated(homonyms)])

#Standardizing municipalities names - lists and data sets
list_mun_std <- correct.mun(list_mun)
mun_concla$Nome_Município <- correct.mun(mun_concla$Nome_Município)
coman_noState$municipality_gbif_std <- correct.mun(coman_noState$municipality_gbif)
coman_state$municipality_gbif_std <- correct.mun(coman_state$municipality_gbif)
coman_noState$county_std <- correct.mun(coman_noState$county)
coman_state$county_std <- correct.mun(coman_state$county)

#Selecting mun_concla municipalities in which the campos rupestres occur
mun_concla_cr <- mun_concla[mun_concla$Nome_Município %in% list_mun_std, ]

#Manually checking and correcting 
#check <- data.frame(plyr::count(coman_noState$municipality_gbif_std)) #ok

#check <- data.frame(plyr::count(coman_noState$county_std)) #ok
coman_na <- rbind(coman_na, 
                       coman_noState[which(coman_noState$municipality_gbif_std %in% 
                                                c("No Disponible")), 1:18])

#check <- data.frame(plyr::count(coman_state$municipality_gbif_std)) #ok
coman_state$municipality_gbif_std[coman_state$municipality_gbif_std == 
                                     "Alto Paraiso"] <- "Alto Paraiso de Goias"
coman_state$municipality_gbif_std[coman_state$municipality_gbif_std == 
                                    "Sao Tome das Letras"] <- "Sao Thome das Letras"
coman_state$municipality_gbif_std[coman_state$municipality_gbif_std == 
                         "Sao Joao Del Rey"] <- "Sao Joao Del Rei"
coman_na <- rbind(coman_na, 
                  coman_noState[which(coman_noState$municipality_gbif_std %in% 
                                        c("Bahia", "Chapada dos Guimaraes",
                                          "Minas Gerais", "Serra do Cipo", "Serra do Cabral")), 1:18])

check <- data.frame(plyr::count(coman_state$county_std)) #ok
coman_state$county_std[coman_state$county_std == 
                                    "Sao Joao Del Rey"] <- "Sao Joao Del Rei"
coman_state$county_std[coman_state$county_std == 
                                    "Sao Tome das Letras"] <- "Sao Thome das Letras"
coman_na <- rbind(coman_na, 
                  coman_state[which(coman_state$county_std %in% 
                                        c("Chapada dos Veadeiros",
                                          "Goias", "No Disponible")), 1:18])


#Filtering coman_noState and coman_state with mun_concla_cr (54 and 649)
coman_noState_filt <- coman_noState %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                                    county_std %in% mun_concla_cr$Nome_Município)
coman_state_filt <- coman_state %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                                county_std %in% mun_concla_cr$Nome_Município)

#=====================================================================================================#

#================================#
# INFERRING MUNICIPALITIES NAMES #
#================================#

#Standardizing 'locality'
coman_na$locality_std <- correct.mun(coman_na$locality)

#Vector with municipalities names from mun_concla to be used with grepl
grepl_munc_concla_cr <- c()
for(i in 1:nrow(mun_concla_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", mun_concla_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(coman_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(coman_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = coman_na$locality_std[j])){
      vec[j] <- mun_concla_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(coman_na$locality_std)] <- NA

#New column with inferred municipality names
coman_na$municipality <- vec

#Removing registers with NA for 'municipality' (122)
coman_na <- coman_na %>% filter(!is.na(municipality))

#=====================================================================================================#

#Concatenating information on municipality into a unique column 
coman_noState_filt$municipality <- NA
for(i in 1:nrow(coman_noState)){
  if(!is.na(coman_noState_filt$municipality_gbif_std[i])){
    coman_noState_filt$municipality[i] <- coman_noState_filt$municipality_gbif_std[i] 
  } else if(!is.na(coman_noState_filt$county_std[i])){
    coman_noState_filt$municipality[i] <- coman_noState_filt$county_std[i]
  }
}

coman_state_filt$municipality <- NA
for(i in 1:nrow(coman_state)){
  if(!is.na(coman_state_filt$municipality_gbif_std[i])){
    coman_state_filt$municipality[i] <- coman_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(coman_state_filt$county_std[i])){
    coman_state_filt$municipality[i] <- coman_state_filt$county_std[i]
  }
}

#Concatenating data sets (825)
coman_noState_filt <- coman_noState_filt[ , -which(names(coman_noState_filt) 
                                                     %in% c("county","county_std",
                                                            "municipality_gbif",
                                                            "municipality_gbif_std"))]

coman_state_filt <- coman_state_filt[ , -which(names(coman_state_filt) 
                                                 %in% c("county","county_std",
                                                        "municipality_gbif",
                                                        "municipality_gbif_std"))]
coman_na <- coman_na[ , -which(names(coman_na) 
                                 %in% c("county", "municipality_gbif", "locality_std"))]

coman_noCoord_inf <- rbind(coman_noState_filt, coman_state_filt, coman_na)

#Registers occurring in homonyms municipalities 
#binding them after removing those without state/province information (824)
reg_hom <- coman_noCoord_inf[coman_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
coman_noCoord_inf <- coman_noCoord_inf %>% filter(!municipality %in% homonyms)
coman_noCoord_inf <- rbind(coman_noCoord_inf, reg_hom) 

rm(coman_na, coman_noState, coman_noState_filt, coman_state, coman_state_filt,
   mun_concla, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, coman_noCoord, reg_hom)

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

#Some municipalities are not included in the data set for centroids
#not_included <- mun_concla_cr %>% filter(!Nome_Município %in% mun_centr_cr$`minor area`)

#Inferring coordinates based on 'municipality'
for(i in 1:length(coman_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(coman_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      coman_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      coman_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Homonyms
for(i in 1:length(coman_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(coman_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       coman_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(coman_noCoord_inf$stateprovince[i])){
      coman_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      coman_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing NA's for coordinates (804)
coman_noCoord_inf <- coman_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, mun_concla_cr, centr_hom, homonyms)

#=====================================================================================================#

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unnimportant columns (2,071)
coman_noCoord_inf <- coman_noCoord_inf %>% dplyr::select(id,
                                                    institutioncode,
                                                    collectioncode,
                                                    catalognumber,
                                                    gen_sp,
                                                    subspecies,
                                                    identifiedby,
                                                    latitude,
                                                    longitude)

coman_coordClean <- coman_coordClean %>% dplyr::select(id,
                                                  institutioncode,
                                                  collectioncode,
                                                  catalognumber,
                                                  gen_sp,
                                                  subspecies,
                                                  identifiedby,
                                                  latitude,
                                                  longitude)

coman_clean <- rbind(coman_coordClean, coman_noCoord_inf)

rm(coman, coman_coordClean, coman_noCoord_inf, i, j, correct.mun, generate.names, replace.names,
   titling)

#Writing *.csv 
#write.csv(coman_clean, file = "Datasets/Comanthera/coman_clean.csv", row.names = FALSE)

#=====================================================================================================#
#====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Reading clean dataset
coman_clean <- read.csv(file = "datasets/Comanthera/coman_clean.csv", na.strings = c("", NA))

#Standardizing gen_sp column
coman_clean$gen_sp <- gsub(" ", "_", coman_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the cr shapefile and, then, with the grids 
coords <- coman_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- coman_clean
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
write.csv(coords, "datasets/Comanthera/coman_cr.csv", row.names = FALSE)
