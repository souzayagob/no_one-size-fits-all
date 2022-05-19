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

#Reading GBIF data
mimosa_gbif <- fread(file = "datasets/Mimosa/0012888-190415153152247_gbif_mimosa/occurrence.txt",
                         na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Reducing data dimensionality by selecting only necessary columns
mimosa_gbif <- mimosa_gbif %>% dplyr::select(institutionCode,
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
mimosa_gbif <- mimosa_gbif %>% rename("species" = specificEpithet,
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
mimosa_gbif <- cbind(id = 1:nrow(mimosa_gbif), mimosa_gbif)

#=============#
# speciesLink #
#=============#

#Reading spLink (Warning message because of double quotes. Not a problem in this context)
mimosa_spLink <- fread(file = "Datasets/Mimosa/speciesLink_mimosa_72660_20190514134517.txt", 
                           na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Selecting important attributes
mimosa_spLink <- mimosa_spLink %>% dplyr::select(institutioncode,
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
mimosa_spLink$longitude <- as.numeric(as.character(mimosa_spLink$longitude))
mimosa_spLink$latitude <- as.numeric(as.character(mimosa_spLink$latitude))

#Giving an unique ID number for each record
mimosa_spLink <- cbind(id = (nrow(mimosa_gbif) + 1):(nrow(mimosa_gbif) + nrow(mimosa_spLink)), mimosa_spLink)

#======================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (116,987), and adding a column to define the original dataset
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

mimosa <- merge.with.source(x = mimosa_gbif,
                                y = mimosa_spLink,
                                name.x = "gbif",
                                name.y = "splink")

rm(merge.with.source, mimosa_gbif, mimosa_spLink)

#======================================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (observations with NA for any of those attributes
#were not considered) (89,918)
a <- mimosa
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
mimosa <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
mimosa <- mimosa %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
mimosa$latitude[mimosa$latitude == 0] <- NA
mimosa$longitude[mimosa$longitude == 0] <- NA

#Removing registers without identifier name (58,930)
mimosa$identifiedby[mimosa$identifiedby %in% c("?", "-", "#", "##", "//05/1987",
                                               "24/X/2013") | mimosa$identifiedby == 0] <- NA
mimosa <- mimosa %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- mimosa$stateprovince
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
mimosa$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
mimosa$stateprovince[mimosa$stateprovince == "?" | mimosa$stateprovince == "-"] <- NA
plyr::count(mimosa$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur 
#according to Silveira et al (2016) (31,575)
mimosa <- mimosa %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                        "Minas Gerais", 
                                                                        "Bahia", 
                                                                        "Pernambuco", 
                                                                        "Paraiba", 
                                                                        "Mato Grosso",
                                                                        "Distrito Federal"))

#Removing records without species level identification (30,750)
mimosa <- mimosa %>% filter(!is.na(species))
mimosa <- mimosa %>% filter(!species %in% c("sp.", "sp1"))
#plyr::count(mimosa$species)

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (30,617)
#plyr::count(mimosa$basisofrecord)
mimosa <- mimosa %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                                 "s", "S"))
mimosa <- mimosa %>% dplyr::select(-basisofrecord)

#=====================================================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

#Creating a vector to work with
identifiedby <- mimosa$identifiedby

#Generating a list of identifiers ordered by frequency of identifications (more identifications
#is supposed to be related with a better identification)
#Counting check
#names.count <- as.data.frame(plyr::count(identifiedby))
#names.count <- names.count[order(-names.count$freq), ]

#R. C. Barneby
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
              check.by = c("R. C. Barneby", "R. C. BARNEBY",
                           "Barneby, R. C", "BARNEBY, R. C."),
              replace.by = "R. C. Barneby",
              not.replace = c("R. Bameby", "R. M. Harley"),
              replace = c("R. Barneby (NY)", "Rupert Barneby",
                          "BArneby", "R. BARNEBY/G. HATSCHBACH",
                          "R. Barneby 85", "R. Barneby (!RL 2014)",
                          "R. Barneby, (!RL 2014)", "R. Barneby (!RL2014)",
                          "R. Barneby 84", "R. Barneby 93", "Rupert C. Barneby",
                          "Rupert Charles Barneby", "R Brandeby", "R. Barneby 1983",
                          "R. Barneby (NY) 1983", "R. Barneby 1986-1991", 
                          "Rupert Charles Barneby; b.1911; d.2000; Barneby", "R. Bameby"))

#L. M. Borges
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                      check.by = c("L. M. Borges", "Borges, L. M.",
                                   "L. M. BORGES", "BORGES, L. M."),
                      replace.by = "L. M. Borges",
              replace = c("Leonardo M. Borges", "L.M. Borges et. al.;",
                          "L.M. Borges (SPF)", "L.M.Borges (SPF)", 
                          "L.M. Borges/13-05-2014"),
              not.replace = c("L. M. G. Nogueira"))

#M. F. Simon
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                      check.by = c("M. F. Simon", "Simon, M. F.",
                                   "M. F. SIMON", "SIMON, M. F."),
                      replace.by = "M. F. Simon",
              not.replace = c("M.S. Simo"),
              replace = c("Marcelo F. Simon", "M. Simon, IN LIT. L. P. Queiroz",
                          "M. Simon & L. M. Borges", "Simon, MF; Queiroz, LP",
                          "Marcelo F Simon", "Marcelo Simon","Simon, MF; Barneby, RC"),
              return = TRUE)

#L. P. de Queiroz
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                      check.by = c("L. P. de Queiroz", "Queiroz, L. P.",
                                   "L. P. DE QUEIROZ", "QUEIROZ, L. P."),
                      replace.by = "L. P. de Queiroz",
              replace = c("Luciano Paganucci de Queiroz", "L.P.de Queiroz & E.R.de Souza",
                          "L.P.de Queiroz & D.Cardoso", "L. P. Queiroz (HUEFS) 2002-01",
                          "L.P.de Queiroz & R.M.Santos", "L.P.de Queiroz & T.S.Nunes",
                          "D. Queiroz", "L. P. Queiroz (HUEFS) 2002-11",
                          "Santos, R.M. dos; Queiroz, L.P. de", "Santos, RM dos; Queiróz, LP de",
                          "Santos, RM; Queiroz, LP de","L. Paganucci", "Queiroz, LP de; Lima, MPM de"),
              not.replace = c("Queiroz, R.T.", "Queiroz, RT"),
                      return = T)

#R. T. de Queiroz
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = c("R. T. de Queiroz", "Queiroz, R. T",
                                           "R. T. DE QUEIROZ", "QUEIROZ, R. T."),
                              replace.by = "R. T. de Queiroz",
                              not.replace = c("L. P. de Queiroz"),
                              replace = c("R.T, Queiroz"),
                              return = T)

#V. F. Dutra
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                      check.by = c("V. F. Dutra", "V. F. DUTRA",
                                   "Dutra, V. F.", "DUTRA, V. F."),
                      replace.by = "V. F. Dutra",
                      replace = c("V.F.Dutra (VIC)", "V.F.Dutra (VIC), (!RL 2014)",
                                  "Valquíria Ferreira Dutra", "Valquiria Ferreira Dutra"),
                      return = T)

#J. S. Silva
identifiedby <- replace.names(x = identifiedby, top = 0.9, bottom = 0.7, 
              check.by = c("J. S. Silva", "Silva, J. S",
                           "J. S. SILVA", "SILVA, J. S."),
              replace.by = "J. S. Silva",
              not.replace = c("Silva, J.P.", "SILVA, R.R.", "Silva, MS",
                              "Silva, J.L.", "Silva, P.", "L. S. Silva",
                              "Silva, J.B.", "Silva, R.", "Silva, MJ" ),
              replace = c("J.S. Silva (2)", "Silva, JS; Nascimento, JGA",
                          "J. S. Silva 2009-08-14", "J. Santos Silva (UEC) 2010-10-27",
                          "J.S. Silva; M. Sales", "Silva, J.S.; Sales, M.", "Santos, J.S.",
                          "J.S. Silva (UEC)", "J. Santos S.", "J. Santos Silva", "J. Santos S. (!RL 2014)",
                          "Juliana Santos Silva (UEC)", "J.S.Silva/14-VIII-2009",
                          "Juliana Santos Silva", "Juliana Santos Silva; Universidade Estadual de Campinas"),
              return = T)

#G. P.  Lewis
identifiedby <- replace.names(x = identifiedby, top = 0.9, bottom = 0.7, 
                      check.by = c("G. P. Lewis", "Lewis, G. P.",
                                   "G. P. LEWIS", "LEWIS, G. P."),
                      replace.by = "G.  P. Lewis",
                      replace = c("Lewis, G.P.; Page, J.S.", "G.P. LEWIS & J.S. PAGE",
                                  "G,P LEWIS", "Lewis, G.P.; Nascimento, M.S.B.",
                                  "Lewis, G.P.; Clark, R.", "G.P. Lewis 78",
                                  "Lewis, G.; Rico, L.", "Gwilym P. Lewis",
                                  "G. P. Lewis & J. S. Page", "E.P. LEWIS", "G.P. Lewis e",
                                  "Lewis, G.P.;Klitgaard, B.", "Lewis, G.P.; Ratter, J.A.",
                                  "G.P. Lewis & J.S. Page", "G. P. Lewis; J. S. Page",
                                  "Lewis, GP; Ratter, JA"))

#M. L. Guedes
identifiedby <- replace.names(x = identifiedby, top = 0.80, bottom = 0.7, 
                      check.by = c("M. L. Guedes", "Guedes, M. L.",
                                   "M. L. GUEDES", "GUEDES, M. L."),
                      replace.by = "M. L. Guedes",
              not.replace = c("M. L. Fonseca"),
              return = TRUE)

#A. P. Savassi-Coutinho
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                      check.by = c("A. P. Savassi-Coutinho", "Savassi-Coutinho, A. P.",
                                   "A. P. SAVASSI-COUTINHO", "SAVASSI-COUTINHO, A.P."),
                      replace.by = "A. P. Savassi-Coutinho",
                      replace = c("Coutinho, A.P.S.", "Coutinho, APS"),
                      return = T)

#L. S. B. Jordão
identifiedby <- replace.names(x = identifiedby, top = 0.9, bottom = 0.7, 
                      check.by = c("L. S. B. Jordão", "Jordão, L. S. B.",
                                   "L. S. B. JORDÃO", "JORDÃO, L. S. B."),
                      replace.by = "L. S. B. Jordão",
                      replace  = c("Lucas S.B. Jordão", "LUCAS S.B.JORDÃO",
                                   "Lucas S. B. Jordão", "LUCAS JORDÃO", "Lucas Josdão",
                                   "Lucas, S.B. Jordão", "L.S.B, Jordão", "Lucas Jordão",
                                   "LucasS.B. Jordão", "Jordao, Lucas", "L.S.B., Jordão",
                                   "L.S. B. Jordão", "L. Jordão","Jordão, L.","L.S.B.Jordão", 
                                   "Jordão, LSB", "L.S.B. Jordão"),
                      return = T)

#O. S. Ribas
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                      check.by = c("O. S. Ribas", "Ribas, O. S.",
                                   "O. S. RIBAS", "RIBAS, O. S."),
                      replace.by = "O. S. Ribas",
                      replace = c("Ribas, OS; Cordeiro, J", "\"\"Ribas, O.S.; Cordeiro, J.\"\"",
                                  "Ribas, OS; Barbosa, E", "O.S. Ribas & J. Codeiro", 
                                  "O.S. Ribas & J. Cordeiro", "Cordeiro, J; Ribas, OS",
                                  "Ribas, OS; Larocca, P", "O.S. Ribas; J. Cordeiro"),
                      return = T)

#I. B. Lima
identifiedby <- replace.names(x = identifiedby, top = 0.93, bottom = 0.7, 
              check.by = c("I. B. Lima", "Lima, I. B.",
                           "I. B. LIMA", "LIMA, I. B."),
              replace.by = "I. B. Lima",
              return = T)

#M. Morales
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                      check.by = c("M. Morales", "Morales, M.",
                                   "M. MORALES", "MORALES, M."),
                      replace.by = "M. Morales",
                      replace = c("Matias Morales", "M. Morols", "m. Morols"),
                      not.replace = c("M.Moraes"),
                      return = T)

#C. W. Fagg
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                      check.by = c("C. W. Fagg", "Fagg, C. W.",
                                   "C. W. FAGG", "FAGG, C. W."),
                      replace.by = "C. W. Fagg",
                      replace = c("Fagg, CW; Borges, LM"),
                      return = T)

#J. G. Nascimento
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                      check.by = c("J. G. Nascimento", "Nascimento, J. G.",
                                   "J. G. NASCIMENTO", "NASCIMENTO, J. G."),
                      replace.by = "J. G. Nascimento",
                      not.replace = c("L. M. Nascimento", "Nascimento, AFS",
                                      "Nascimento, FHF", "Nascimento, L.M.",
                                      "Nascimento, M.S.B.", "J.M. Nascimento",
                                      "Nascimento"),
                      replace = c("J. G. A. do Nascimento", "J.G.A.do Nascimento",
                                  "J.G.A. do Nascimento"),
                      return = T)

#R. Grether
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                      check.by = c("R. Grether", "Grether, R.",
                                   "R. GRETHER", "GRETHER, R."),
                      replace.by = "R. Grether",
                      return = T)

#E. Córdula
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = c("E. Córdula", "Córdula, E.",
                                           "E. CÓRDULA", "CÓRDULA, E"),
                              replace.by = "E. Córdula",
                              return = T)

#A. Bocage
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = c("A. Bocage", "Bocage, A.",
                                           "A. BOCAGE", "BOCAGE, A."),
                              replace.by = "A. Bocage",
                              replace = c("A.L.du Bocage Lima", "Bocage, A., Marques, J.S."),
                              return = T)

#J. R. Pirani
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = c("J. R. Pirani", "Pirani, J. R.",
                                           "J. R. PIRANI", "PIRANI, J. R."),
                              replace.by = "J. R. Pirani",
                              replace = c("Pirani, JR; Siniscalchi, CM"),
                              return = T)

#R. H. Fortunato
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = generate.names("R. H.", "Fortunato"),
                              replace.by = "J. R. Pirani",
                              replace = c("Reenée Fortunato", "Renée H. Fortunato",
                                          "Renee H. Fortunato", "Renée H. Forunato"),
                              return = T)

#H. C. de Lima
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = generate.names("H. C.", "de Lima"),
                              replace.by = "H. C. de Lima",
                              not.replace = c("R .C. de Lima"),
                              replace = c("H.C. Lima", "Lima, HC",
                                          "Lima, H.C.", "H. C. de Lima & L.F.G da Silva",
                                          "H C Lima", "H.C. Lima & C.M.J. Mattos",
                                          "H.C.LIMA", "H.C. Lima , C. M. B. Correia",
                                          "H.C. de LIma", "H.C. de Lima & Marli Pires",
                                          "Haroldo C. de Lima", "H.C.Lima", "H. C. de Lima & L. F. G. da Silva"),
                              return = T)

#A. Burkart
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = generate.names("A.", "Burkart"),
                              replace.by = "A. Burkart",
                              replace = c("A. Burkart (!RL 2014)", "Buskart",
                                          "A. E. Burkart; L. B. Smith"),
                              return = T)

#A. M. Miranda
identifiedby <- replace.names(x = identifiedby, top = 0.90, bottom = 0.7, 
                              check.by = generate.names("A. M.", "Miranda"),
                              replace.by = "A. M. Miranda",
                              replace = c("A.M.Miranda et H.Gomes", "A.M. Miranda; M.L. Guedes"),
                              return = T)

#Counting check
names.count <- as.data.frame(plyr::count(identifiedby))
names.count <- names.count[order(-names.count$freq), ]

#Replacing column
mimosa$identifiedby <- identifiedby

#Filtering (20,672)
specialists <- c("R. C. Barneby",
                 "L. M. Borges",
                 "L. P. de Queiroz",
                 "R. T. de Queiroz",
                 "J. S. Silva",
                 "M. S. Simon",
                 "G. P. Lewis",
                 "L. S. B. Jordão",
                 "A. P. Savassi-Coutinho",
                 "M. L. Guedes",
                 "O. S. Ribas",
                 "J. G. Nascimento",
                 "J. R. Pirani",
                 "I. B. Lima",
                 "M. Morales",
                 "C. W. Fagg",
                 "R. Grether",
                 "E. Córdula",
                 "A. Bocage",
                 "R. H. Fortunato",
                 "H. C. de Lima",
                 "A. Burkart",
                 "A. M. Miranda")
mimosa <- mimosa %>% filter(identifiedby %in% specialists)

rm(identifiedby, specialists)
#Creating *.csv
#write.csv(mimosa, file = "Datasets/Mimosa/mimosa.csv", row.names = FALSE)
#======================================================================================#
#Loading *.csv
#mimosa <- fread(file = "Datasets/Mimosa/mimosa.csv",na.strings = c("", NA))

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

#Generating a column for scientific names (without authors and including infraspecific epithet)
mimosa$gen_sp <- paste(mimosa$genus,
                       mimosa$species,
                       mimosa$subspecies,
                       sep = " ")

#taxa <- plyr::count(mimosa$gen_sp)
#taxa$x <- as.character(taxa$x)
#taxa_suggested <- tibble("taxa" = taxa$x, "suggested" = NA)

#Removing NA (character derived from 'subspecies' attribute)
#for(i in 1:nrow(taxa_suggested)){
 # taxa_suggested[i, 1] <- gsub("NA", "", taxa_suggested[i, 1])
  #taxa_suggested[i, 1] <- trimws(taxa_suggested[i, 1])
#}

#Suggesting with flora
#for(i in 1:nrow(taxa_suggested)){
#taxa_suggested$suggested[i] <- suggest.names(taxa_suggested$taxa[i])
#}

#Writing *.csv for manual checking
#write.csv(taxa_suggested, file = "lists/Mimosa/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Mimosa/taxa_corrected.csv", stringsAsFactors = F, 
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

#Removing invalid taxa (20,628)
mimosa$gen_sp <- gsub("NA", "", mimosa$gen_sp)
mimosa$gen_sp <- trimws(mimosa$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen) | taxa_gensp$gen != "Mimosa"]
mimosa <- mimosa[!mimosa$gen_sp %in% invalid_taxa, ]

#Correcting the dataset
taxa_gensp$replace <- taxa_corrected$taxa
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(mimosa)){
  for(j in 1:nrow(taxa_gensp)){
    if(mimosa$gen_sp[i] == taxa_gensp$replace[j]){
      mimosa$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

#Last corrections
#mimosa$gen_sp[!mimosa$gen_sp %in% paste(taxa_gensp$gen, taxa_gensp$sp, sep = " ")]
mimosa$gen_sp[mimosa$gen_sp == "Mimosa ourobrancoã«nsis"] <- "Mimosa ourobrancoensis"

rm(taxa_corrected, taxa_gensp, invalid_taxa, str)

#Creating *.csv
#write.csv(mimosa, file = "Datasets/Mimosa/mimosa.csv", row.names = FALSE)
#======================================================================================#
#Loading *.csv
#mimosa <- fread(file = "Datasets/Mimosa/mimosa.csv",na.strings = c("", NA))

#================#
# SPLITTING DATA #
#================#

#Coords (9,768)
mimosa_coord <- mimosa %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (10,860)
mimosa_noCoord <- mimosa %>% filter(is.na(latitude) | is.na(longitude))

#Two registers, ids 115454 e 115984, have invalid coordinates
#Removing them and transfering to mimosa_noCoord
mimosa_coord <- mimosa_coord %>% filter(!id %in% c(115984, 115454))
mimosa_noCoord <- rbind(mimosa_noCoord, mimosa[mimosa$id %in% c(115984, 115454), ])

#Writing *.csv
#write.csv(mimosa_coord, file = "Datasets/Mimosa/mimosa_coord.csv", row.names = FALSE)
#write.csv(mimosa_noCoord, file = "Datasets/Mimosa/mimosa_noCoord.csv", row.names = FALSE)
#======================================================================================#
#Loading *.csv
#mimosa_coord <- fread(file = "Datasets/Mimosa/mimosa_coord.csv", na.strings = c("", NA))

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (9,546)
mimosa_coordFlagged <- mimosa_coord %>% clean_coordinates(lon = "longitude",
                                                        lat = "latitude",
                                                        species = "gen_sp",
                                                        value = "flagged",
                                                        tests = c("equal", "gbif", 
                                                                  "institutions", 
                                                                  "outliers", "seas",
                                                                  "zeros"))
                                                                  

invalid_coords <- mimosa_coord[mimosa_coordFlagged == FALSE, ]
mimosa_coordClean <- mimosa_coord[mimosa_coordFlagged  == TRUE, ]

#Binding invalid_coords to mimosa_noCoord (11,082)
mimosa_noCoord <- rbind(mimosa_noCoord, invalid_coords)

rm(invalid_coords, mimosa_coordFlagged, mimosa_coord)

#Writing *.csv
#write.csv(mimosa_coordClean, file = "Datasets/Mimosa/mimosa_coordClean.csv", row.names = FALSE)
#write.csv(mimosa_noCoord, file = "Datasets/Mimosa/mimosa_noCoord.csv", row.names = FALSE)
#======================================================================================#
#Loading *.csv
#mimosa_noCoord <- fread(file = "Datasets/Mimosa/mimosa_noCoord.csv", na.strings = c("", NA))

#===============================#
# CLEANING MUNICIPALITIES NAMES #
#===============================#

#Dividing dataset in registers with and without values for stateprovince 
mimosa_noState <- mimosa_noCoord %>% filter(is.na(stateprovince))
mimosa_state <- mimosa_noCoord %>% filter(!is.na(stateprovince))

#NA values for municipalities 
mimosa_na <- mimosa_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(mimosa_noState$county))
#names.count <- names.count[order(-names.count$freq), ]

#Reading list of municipalities according to DTB (Divisão Territorial Brasileira - IBGE) - ftp://geoftp.ibge.gov.br/organizacao_do_territorio/estrutura_territorial/divisao_territorial/
mun_concla <- fread(file = "lists/RELATORIO_DTB_BRASIL_MUNICIPIO.csv",
                    na.strings = c("", NA), encoding = "UTF-8")
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
mimosa_noState$municipality_gbif_std <- correct.mun(mimosa_noState$municipality_gbif)
mimosa_state$municipality_gbif_std <- correct.mun(mimosa_state$municipality_gbif)
mimosa_noState$county_std <- correct.mun(mimosa_noState$county)
mimosa_state$county_std <- correct.mun(mimosa_state$county)

#Selecting mun_concla municipalities in which the campos rupestres occur
mun_concla_cr <- mun_concla[mun_concla$Nome_Município %in% list_mun_std, ]

#nrow(mun_concla_cr) is not equal to length(list_mun_std). Why?
#Municipalities from the São Paulo state
#list_mun_std[!list_mun_std %in% mun_concla_cr$Nome_Município]

#Manually checking and correcting
#check <- data.frame(plyr::count(mimosa_noState$municipality_gbif_std)) #ok
mimosa_noState$municipality_gbif_std[mimosa_noState$municipality_gbif_std == 
                                     "Gouvea"] <- "Gouveia"
mimosa_noState$municipality_gbif_std[mimosa_noState$municipality_gbif_std == 
                                     "Anaje"] <- "Anage"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Sao Joao D'Alianca"] <- "Sao Joao D'alianca"
mimosa_na <- rbind(mimosa_na, 
                   mimosa_noState[which(mimosa_noState$municipality_gbif_std %in% 
                                        c("Brasil", "Goias")), 1:18])

#check <- data.frame(plyr::count(mimosa_noState$county_std)) #ok
mimosa_na <- rbind(mimosa_na, 
                   mimosa_noState[which(mimosa_noState$county_std %in% 
                                          c("Brasil")), 1:18])

#check <- data.frame(plyr::count(mimosa_state$municipality_gbif_std))
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Alto Paraiso"] <- "Alto Paraiso de Goias"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Anaje"] <- "Anage"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Barra de Mendes"] <- "Barra do Mendes"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Betania/Floresta"] <- "Betania"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Corumba  de Goias"] <- "Corumba de Goias"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Gouvea"] <- "Gouveia"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Planaltina de Goias"] <- "Planaltina"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Rio das Contas"] <- "Rio de Contas"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Sao Tome das Letras"] <- "Sao Thome das Letras"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Sao Joao da Alianca"] <- "Sao Joao D'alianca"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std 
                                     %in% c("Vila Bela Santissima Trindade", "Vila Bela da Sma Trindade",
                                            "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Catugi"] <- "Catuji"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Cajazeiras Sao Jose das Piranhas"] <- "Cajazeiras"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Felixandia"] <- "Felixlandia"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Santana do Pirapama"] <- "Santana de Pirapama"
mimosa_state$municipality_gbif_std[mimosa_state$municipality_gbif_std == 
                                     "Terezina de Goias"] <- "Teresina de Goias"
mimosa_na <- rbind(mimosa_na, 
                   mimosa_state[which(mimosa_state$municipality_gbif_std %in% 
                                        c("Chapada da Cotagem", "Chapada Diamantina",
                                          "Chapada do Araripe", "Chapada dos Guimaraes",
                                          "Chapada dos Veadeiros", "Paraiba", "Pernambuco",
                                          "Serra do Cipo", "Serra do Espinhaco", "Serra do Tombador",
                                          "Ufba Ondina", "Veadeiros", "Bahia", "Chapadao do Ceu", "Goias", 
                                          "Minas Gerais", "Serra do Cipo", "Serra do Espinhaco", 
                                          "Serra do Tombador", "Veadeiros")), 1:18])

#check <- data.frame(plyr::count(mimosa_state$county_std))
mimosa_na <- rbind(mimosa_na, 
                   mimosa_state[which(mimosa_state$county_std %in% 
                                        c("", "Chapada dos Guimaraes",
                                          "Chapada dos Veadeiros", "Chapada Gaucha",
                                          "Goias", "Minas Gerais", "No Disponible",
                                          "Pe", "Ufba Ondina", "Veadeiros", 
                                          "Minas Gerais", "Mato Grosso", "Salvador")), 1:18])
mimosa_state$county_std[mimosa_state$county_std %in% 
                                     c("Alto Paraiso", "Alto Paraa-so de Goias")] <- "Alto Paraiso de Goias"
mimosa_state$county_std[mimosa_state$county_std == 
                                     "Anaje"] <- "Anage"
mimosa_state$county_std[mimosa_state$county_std == 
                                     "Barra de Mendes"] <- "Barra do Mendes"
mimosa_state$county_std[mimosa_state$county_std == 
                                     "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
mimosa_state$county_std[mimosa_state$county_std == 
                          "Cristalina Mun"] <- "Cristalina"
mimosa_state$county_std[mimosa_state$county_std == 
                          "Joaquim Fela-cio"] <- "Joaquim Felicio"
mimosa_state$county_std[mimosa_state$county_std == 
                          "Rio das Contas"] <- "Rio de Contas"
mimosa_state$county_std[mimosa_state$county_std %in% 
                          c("Sao Joao da Alianca", 
                            "Sao Joao D##Alianca")] <- "Sao Joao D'alianca"
mimosa_state$county_std[mimosa_state$county_std %in% 
                          c("Vila Bela de Santa-ssima Trindade", 
                            "Vila Bela de Santissima Trindade",
                            "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
mimosa_state$county_std[mimosa_state$county_std == 
                          "Santana de Pirapama"] <- "Santana do Pirapama"
mimosa_state$county_std[mimosa_state$county_std == 
                          "Brasa-lia"] <- "Brasilia"

#Filtering mimosa_noState and mimosa_state with mun_concla_cr
mimosa_noState_filt <- mimosa_noState %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                                    county_std %in% mun_concla_cr$Nome_Município)
mimosa_state_filt <- mimosa_state %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                                county_std %in% mun_concla_cr$Nome_Município)

#================================#
# INFERRING MUNICIPALITIES NAMES #
#================================#

#Standardizing 'locality'
mimosa_na$locality_std <- correct.mun(mimosa_na$locality)

#Vector with municipalities names from mun_concla to be used with grepl
grepl_munc_concla_cr <- c()
for(i in 1:nrow(mun_concla_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", mun_concla_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(mimosa_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(mimosa_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = mimosa_na$locality_std[j])){
      vec[j] <- mun_concla_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(mimosa_na$locality_std)] <- NA

#New column with inferred municipality names
mimosa_na$municipality <- vec

#Removing registers with NA for 'municipality'
mimosa_na <- mimosa_na %>% filter(!is.na(municipality))

#Concatenating information on municipality into a unique column 
mimosa_noState_filt$municipality <- NA
for(i in 1:nrow(mimosa_noState)){
  if(!is.na(mimosa_noState_filt$municipality_gbif_std[i])){
    mimosa_noState_filt$municipality[i] <- mimosa_noState_filt$municipality_gbif_std[i] 
  } else if(!is.na(mimosa_noState_filt$county_std[i])){
    mimosa_noState_filt$municipality[i] <- mimosa_noState_filt$county_std[i]
  }
}

mimosa_state_filt$municipality <- NA
for(i in 1:nrow(mimosa_state)){
  if(!is.na(mimosa_state_filt$municipality_gbif_std[i])){
    mimosa_state_filt$municipality[i] <- mimosa_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(mimosa_state_filt$county_std[i])){
    mimosa_state_filt$municipality[i] <- mimosa_state_filt$county_std[i]
  }
}

#Concatenating data sets
mimosa_noState_filt <- mimosa_noState_filt %>% 
  dplyr::select(colnames(mimosa_noState_filt)[!colnames(mimosa_noState_filt) %in%
                                                    c("county","county_std",
                                                      "municipality_gbif",
                                                      "municipality_gbif_std")])

mimosa_state_filt <- mimosa_state_filt %>% 
  dplyr::select(colnames(mimosa_state_filt)[!colnames(mimosa_state_filt) %in%
                                                  c("county","county_std",
                                                    "municipality_gbif",
                                                    "municipality_gbif_std")])

mimosa_na <- mimosa_na %>% 
  dplyr::select(colnames(mimosa_na)[!colnames(mimosa_na) %in%
                                          c("county","county_std",
                                            "municipality_gbif",
                                            "municipality_gbif_std")])

mimosa_noCoord_inf<- rbind(mimosa_noState_filt, mimosa_state_filt, mimosa_na, fill = TRUE)

#Registers occurring in homonyms municipalities
reg_hom <- mimosa_noCoord_inf[mimosa_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
mimosa_noCoord_inf <- mimosa_noCoord_inf %>% filter(!municipality %in% homonyms)
mimosa_noCoord_inf <- rbind(mimosa_noCoord_inf, reg_hom) #binding them after removing those without state/province information

rm(mimosa_na, mimosa_noState, mimosa_noState_filt, mimosa_state, mimosa_state_filt,
   mun_concla, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, mimosa_noCoord, reg_hom)

#Writing *.csv 
#write.csv(mimosa_noCoord_inf, file = "Datasets/Mimosa/mimosa_noCoord_inf.csv", row.names = FALSE)
#======================================================================================#
#Loading *.csv
#mimosa_noCoord_inf <- fread(file = "Datasets/Mimosa/mimosa_noCoord_inf.csv", na.strings = c("", NA))

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
for(i in 1:length(mimosa_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(mimosa_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      mimosa_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      mimosa_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Homonyms
for(i in 1:length(mimosa_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(mimosa_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       mimosa_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(mimosa_noCoord_inf$stateprovince[i])){
      mimosa_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      mimosa_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing NA's for coordinates
mimosa_noCoord_inf <- mimosa_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, mun_concla_cr, centr_hom, homonyms)

#Writing *.csv 
#write.csv(mimosa_noCoord_inf, file = "Datasets/Mimosa/mimosa_noCoord_inf.csv", row.names = FALSE)
#======================================================================================#
#Loading *.csv
#mimosa_noCoord_inf <- fread(file = "Datasets/Mimosa/mimosa_noCoord_inf.csv", na.strings = c("", NA))
#mimosa_coordClean <- fread(file = "Datasets/Mimosa/mimosa_coordClean.csv", na.strings = c("", NA))

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unimportant columns
mimosa_noCoord_inf <- mimosa_noCoord_inf %>% dplyr::select(id,
                                                    institutioncode,
                                                    collectioncode,
                                                    catalognumber,
                                                    gen_sp,
                                                    subspecies,
                                                    identifiedby,
                                                    latitude,
                                                    longitude)

mimosa_coordClean <- mimosa_coordClean %>% dplyr::select(id,
                                                  institutioncode,
                                                  collectioncode,
                                                  catalognumber,
                                                  gen_sp,
                                                  subspecies,
                                                  identifiedby,
                                                  latitude,
                                                  longitude)

mimosa_clean <- rbind(mimosa_coordClean, mimosa_noCoord_inf)

#Is there any duplicated registers? 
#mimosa_clean[duplicated(mimosa_clean$id) == FALSE, ] #no

rm(mimosa, mimosa_coordClean, mimosa_noCoord_inf, i, j, correct.mun, generate.names,
   replace.names, titling)

#Writing *.csv 
#write.csv(mimosa_clean, file = "Datasets/Mimosa/mimosa_clean.csv", row.names = FALSE)

#====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Reading clean dataset
mimosa_clean <- read.csv(file = "datasets/Mimosa/mimosa_clean.csv", na.strings = c("", NA))

#Standardizing gen_sp column
mimosa_clean$gen_sp <- gsub(" ", "_", mimosa_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the cr shapefile and, then, with the grids 
coords <- mimosa_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- mimosa_clean
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
write.csv(coords, "datasets/Mimosa/mimosa_cr.csv", row.names = FALSE)
