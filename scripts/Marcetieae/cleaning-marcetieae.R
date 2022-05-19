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

#Reading GBIF data (24,300)
marcetieae_gbif <- fread(file = "datasets/Marcetieae/0245858-200613084148143/occurrence.txt",
                         na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Reducing data dimensionality by selecting only necessary columns
marcetieae_gbif <- marcetieae_gbif %>% dplyr::select(institutionCode,
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
marcetieae_gbif <- marcetieae_gbif %>% rename("species" = specificEpithet,
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
marcetieae_gbif <- cbind(id = 1:nrow(marcetieae_gbif), marcetieae_gbif)

#=============#
# speciesLink #
#=============#

#Reading spLink (21,924)
marcetieae_spLink <- fread(file = "datasets/Marcetieae/speciesLink_all_93649_20210408174219.txt", 
                           na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Selecting important attributes
marcetieae_spLink <- marcetieae_spLink %>% dplyr::select(institutioncode,
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
marcetieae_spLink$longitude <- as.numeric(as.character(marcetieae_spLink$longitude))
marcetieae_spLink$latitude <- as.numeric(as.character(marcetieae_spLink$latitude))

#Giving an unique ID number for each record
marcetieae_spLink <- cbind(id = (nrow(marcetieae_gbif) + 1):(nrow(marcetieae_gbif) + nrow(marcetieae_spLink)), marcetieae_spLink)

#=====================================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (46,224), and adding a column to define the original dataset
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

marcetieae <- merge.with.source(x = marcetieae_gbif,
                                y = marcetieae_spLink,
                                name.x = "gbif",
                                name.y = "splink")

rm(merge.with.source, marcetieae_gbif, marcetieae_spLink)

#=====================================================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (observations with NA for any of those attributes
#were not considered) (38,201)
a <- marcetieae
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
marcetieae <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
marcetieae <- marcetieae %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
marcetieae$latitude[marcetieae$latitude == 0] <- NA
marcetieae$longitude[marcetieae$longitude == 0] <- NA

#Removing registers without determiner's name (3,620) 
plyr::count(marcetieae$identifiedby)
marcetieae$identifiedby[marcetieae$identifiedby %in% c("?", "-", "#", "##", "//05/1987",
                                                       "24/X/2013", "`", "15") | is.numeric(marcetieae$identifiedby)] <- NA
marcetieae <- marcetieae %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- marcetieae$stateprovince
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
marcetieae$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
marcetieae$stateprovince[marcetieae$stateprovince == "?" | marcetieae$stateprovince == "-"] <- NA
plyr::count(marcetieae$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur (13,448)
#according to Silveira et al (2016) 
marcetieae <- marcetieae %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                                "Minas Gerais", 
                                                                                "Bahia", 
                                                                                "Pernambuco", 
                                                                                "Paraiba", 
                                                                                "Mato Grosso",
                                                                                "Distrito Federal"))

#Removing records without species level identification (13,049)
marcetieae <- marcetieae %>% filter(!is.na(species))

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (13,041)
plyr::count(marcetieae$basisofrecord)
marcetieae <- marcetieae %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                                         "s", "S", "PreservedSpecimen", "PreserverdSpecimen"))
marcetieae <- marcetieae %>% dplyr::select(-basisofrecord)

#=====================================================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

#Creating a vector to work with
identifiedby <- marcetieae$identifiedby

#Extracting a vector including determiners' names. It is preferable to use this vector instead of the data set itself because any changes are easily reversible. 
identifiedby <- as.character(marcetieae$identifiedby)

#Ideally, it is better to start with a list of taxonomic specialists' compiled beforehand. Alternatively, as experts likely indentified the majority of samples from a given taxon, it is possible to infer specialists based on identification frequency. In this example we looked for specialists at the top of the list below. 
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#To improve accuracy, we confirmed if names in 'names.count' with at least 3 identifications were specialists by searching for taxonomic publications for the family of the focal group and authored by each name at Google Scholar. In this example, we searched for: allintitle: Melastomataceae OR Marcetieae author:"determiner".

#Next, based on the function 'replace.names', we standardized specialist's name. This is done in two iterations:
#(1) The first iteration returns, for manual evaluation, the automatically replaced names (names above the 'top' threshold) and names that are worth to be checked (names above the 'bottom' threshold but below the 'top' threshold).
#(2) In the second iteration, names that were erroneously replaced in the first iteration should be included in the argument 'not replace'. Likewise, names that were supposed to be replaced but were below the 'top' threshold should be included in the argument 'replace'.

#Because the procedure is iterative, the user should not change the vector 'identifiedBy' directly in the first iteration. For this purpose, we used a secondary vector ('identifiedBy_2'). After the second iteration, when everything should be set, the user must assign 'identifiedBy' as 'identifiedBy_2' before running the protocol for the following name. 

#At the end of the standardizing procedure, the following vector should contain all specialists' names
specialists <- c()

#JJ Wurdack
replace.by <- "JJ Wurdack"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("J. J. Wurdack 1957","Wurdack, J. J. 1973, , Rev. Freire-Fierro, A.",
                                            "J. J. Wurdack, 85","J. J. Wurdack/Dez.1990","Wurdack/1993",
                                            "Wurdack/1991","John Julius Wurdack; b.1921; d.1988; Wurdack, J.J."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#R Romero
replace.by <- "R Romero"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("R. Romfro 2002-08","R. romero"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#AB Martins
replace.by <- "AB Martins"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Martins, E."),
                                replace = c("Angela Borges Martins (por foto)",
                                            "A. B. Martins (UEC) 1997-07",
                                            "A. B. Martins (UEC) 1990-04",
                                            "Ângela B. Martins",
                                            "Semir, J; Martins, AB",
                                            "A. B. Martins (UEC) 1990-02",
                                            "A. B. Martins 1996-08"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#PJF Guimarães
replace.by <- "PJF Guimarães"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#AKA Santos
replace.by <- "AKA Santos"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Santos, C.P.","Santos, TM","Santos, R. M.",
                                                "Santos R.M.","Santos, T.S.","Santos, R.M.",
                                                "Santos, E.M.","Santos, EM","Sant'Ana, AK"),
                                replace = c("A. K. A. dos Santos","Andréa Karla Almeida dos Santos",
                                            "A.K.A.Santos & A.B.Martins","A. Karla Santos",
                                            "Andrea K A Santos","A. K. A. santos","Aka Santos",
                                            "J.G.Freitas & A.K.A.Santos"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#A Freire-Fierro
replace.by <- "A Freire-Fierro"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Alina Freire Fierro","Rev. Freire-Fierro, A.",
                                            "Alina Freire-Fierro","Monografia Freire-Fierro",
                                            "Alina Freire-Fierro; PH Herbarium; Botany; b.1964; Freire-Fierro",
                                            "Fierro, A.F."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#D Nunes da Silva
replace.by <- "D Nunes da Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("D. Nunes","Nunes, D","D.NUNES",
                                            "Nunes, D.","D.Nunes"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#SAC Chiea
replace.by <- "SAC Chiea"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("S.A.C. Chiea; I.B.S.P.",
                                            "Chiea, S. C. 10/1995, Rev. Freire-Fierro, A."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JFA Baumgratz
replace.by <- "JFA Baumgratz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Guedes, ML; Baumgratz, JFA",
                                            "Guedes, M. L. | Baumgratz, J. F. A.",
                                            "J. F. A. Baumgratz & M. L. Guedes",
                                            "J.F. Baumgratz & M.L. Guedes",
                                            "J. F. Baumgratz & M. L. Guedes",
                                            "J.F.A.Baumgratz & M.R.Moutinho",
                                            "J.F.A. Baumgratz & B. Valente/12-03-2018",
                                            "J. F. A. Braumgratz","J.F.A. Baumgratz & B. Valente",
                                            "J. F. A. Baumgratz & P. Rosa","J.F.A.BAUMGRATZ & L.P.G.ROSA",
                                            "Guedes, M.L.; Baumgratz, J.F.A.","J.F.A. Baumgratz et D.S.P. Silva",
                                            "J.F.A. Baumgratz & K.C. Silva-Gonçalves","J.F.A. baumgratz & B, Chiavegatto",
                                            "J.F.A. Baumgratz & B. Chiavegatto","J.F.A.Baungratz, M.L.Souza & B.Chiavegatto",
                                            "J.F. BAUMGRATZ","J.F. BAUMGRATZ & M.L. GUEDES","J.F.A.Baumgratz e L.A.F.Santos Filho",
                                            "J.F.A. Baumgratz & B. Chiavegattoi","J. F. A. Baumgratz & K. C. Silva",
                                            "J.F.A. Braumgratz & B. Chiavegatto","REV: Baumgratz, J. F. A. & Santos Filho, L. A. F.",
                                            "Guedes, MLS; Baumgratz, JFA","Guedes, M.L.; Baumgratz, J.F.","Guedes, M.L.; Baumgratz, A.",
                                            "Souza, M.L.; Baumgratz, J.F.A.","Guedes, ML; ML; Baumgratz, JFA; Amaral Júnior, A"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#DN da Silva
replace.by <- "DN da Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("D.R. da Silva"),
                                replace = c("Silva, D.N."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#F Almeda
replace.by <- "F Almeda"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("F. Almeda & O. Robinson","F. Almeda, O. Robinson",
                                            "F. Almeda & O. Robinsn","F. ALmeda & O. Robinson",
                                            "Almeda & O. Robinson 1989","F. Almeda & O Robinson",
                                            "F. Almeda e O. Robinson","F. Almeda, & O. Robinson",
                                            "f. Almeda & O. Robinson","F. Almeda & O. Roninson",
                                            "F. Almeda; O. Robinson","Almeda & Robinson",
                                            "F. Almedo & O. Robinson","F. almeda & O. Robinson",
                                            "F. Almeda & O. robinson","F. Almeda & O. Rinson",
                                            "F. Ameda & O. Robinson","F. Almeda & O. Robinson (CAS)",
                                            "F. Almeda & O. Robinson 1989","Almeda & Robinson A.",
                                            "Robinson, O; Almeda, F"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#SS Renner
replace.by <- "SS Renner"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("L. Renner", "Rennel, J.I."),
                                replace = c("Susanne Sabine Renner; b.1954; S.S.Renner",
                                            "S. S. Renner 1986","NYBG 50:84, Renner Mem.",
                                            "Susanne Sabine Renner"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#R Goldenberg
replace.by <- "R Goldenberg"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Semir, J; Goldenberg, R","Renato Goldenberg","Morais, J.W. & Goldenberg, R."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#EM Woodgyer
replace.by <- "EM Woodgyer"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#CBR Munhoz
replace.by <- "EM Woodgyer"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#CBR Munhoz
replace.by <- "CBR Munhoz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#MJR Rocha
replace.by <- "MJR Rocha"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Rocha, V."),
                                replace = c("Ribeiro, RC; Rocha, MJR","M. J. R. Rocha; G. W. Fernandes",
                                            "M.J.R. Rocha; L.S. Leoni"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#AFA Versiane
replace.by <- "AFA Versiane"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("A. F. A. Versiane; M. J. R. Rocha",
                                            "A.F.A. Versiane; M.J.R.R. Rocha"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#RC Seco
replace.by <- "RC Seco"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#IM Araújo
replace.by <- "IM Araújo"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Araújo, CMLR","Araujo, CMLR","Araújo, C.M.R.L.",
                                                "Araújo, D","Araújo, CMR","Araújo, J.A.G.",
                                                "Araújo, C.M.L.R."),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#CMLR Araújo
replace.by <- "CMLR Araújo"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("IM Araújo","Araújo, D","Araújo, J.A.G."),
                                replace = c("C. M. L. R. Araújo"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#J Coelho de Jesus
replace.by <- "J Coelho de Jesus"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("COELHO,J.","Coelho, J.O.","Coelho, J.",
                                            "J.Coelho","J. Coelho"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#FA Michelangeli
replace.by <- "FA Michelangeli"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("F. A. Michelangeli; R. Goldenberg",
                                            "FABIAN A. MICHELANGELI"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#B Chiavegatto
replace.by <- "B Chiavegatto"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("B.Eliavegato"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#R Kriebel
replace.by <- "R Kriebel"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LL Justino 
replace.by <- "LL Justino"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#FC Hoehne
replace.by <- "FC Hoehne"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#TP Rolim
replace.by <- "TP Rolim"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#SG da Vinha
replace.by <- "SG da Vinha"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#ML Souza
replace.by <- "ML Souza"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Souza, V.C.","Souza, T.","Souza, VC",
                                                "Sousa, T","Souza, J.P.","Souza, I.","Souza, CSD",
                                                "L.F. Souza","Souza, V. C."),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#FS Meyer
replace.by <- "FS Meyer"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Meyer, P.B."),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#KF Rodrigues
replace.by <- "KF Rodrigues"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("K.F. Rodrigues; C.P. Candido",
                                            "Fidanza Rodrigues, K; Martins, A; Rocha, MJR",
                                            "Fidanza Rodrigues, K"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#AR de Rezende
replace.by <- "AR de Rezende"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Rezende, AR","A.R.Rezende","Rezende, A.R."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#AM Amorim
replace.by <- "AM Amorim"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#BDN Valente
replace.by <- "BDN Valente"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Valente, G.E.","Vale, F"),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#Replacing column
marcetieae$identifiedby <- identifiedby

#Filtering (9,940)
marcetieae <- marcetieae %>% filter(identifiedby %in% specialists)

rm(identifiedby, identifiedby_2, specialists, names.count, replace.by)

#=====================================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

library(flora)

#Generating a column for scientific names (without authors and including infraspecific epithet)
marcetieae$gen_sp <- paste(marcetieae$genus,
                           marcetieae$species,
                           marcetieae$subspecies,
                           sep = " ")

#Names
taxa <- plyr::count(marcetieae$gen_sp)
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
write.csv(taxa_suggested, file = "lists/Marcetieae/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Marcetieae/taxa_corrected.csv", stringsAsFactors = F, 
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

#Removing invalid taxa (9,895)
marcetieae$gen_sp <- gsub("NA", "", marcetieae$gen_sp)
marcetieae$gen_sp <- trimws(marcetieae$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen) | !taxa_gensp$gen %in% c("Acanthella",
                                                                                  "Aciotis",
                                                                                  "Acisanthera",
                                                                                  "Dicrananthera",
                                                                                  "Noterophila",
                                                                                  "Rostranthera",
                                                                                  "Appendicularia",
                                                                                  "Comolia",
                                                                                  "Leiostegia",
                                                                                  "Ernestia",
                                                                                  "Fritzschia",
                                                                                  "Macairea",
                                                                                  "Marcetia",
                                                                                  "Nepsera",
                                                                                  "Pseudoernestia",
                                                                                  "Sandemania",
                                                                                  "Siphanthera",
                                                                                  "Comoliopsis",
                                                                                  "Brasilianthus")]
marcetieae <- marcetieae[!marcetieae$gen_sp %in% invalid_taxa, ]

#Correcting the data set (varieties suppressed)
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(marcetieae)){
  for(j in 1:nrow(taxa_gensp)){
    if(marcetieae$gen_sp[i] == taxa_gensp$replace[j]){
      marcetieae$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

rm(taxa, taxa_suggested, taxa_corrected, taxa_gensp, invalid_taxa, str)

#=====================================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (4,510)
marcetieae_coord <- marcetieae %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (5,385)
marcetieae_noCoord <- marcetieae %>% filter(is.na(latitude) | is.na(longitude))

#Removing invalid coordinates and transferring them to marcetieae_noCoord (invalid coordinates are indicated by function clean_coordinates() in the next session)
marcetieae_noCoord <- rbind(marcetieae_noCoord, marcetieae_coord[c(), ])
marcetieae_coord <- marcetieae_coord[-c(), ]

#=====================================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (4,309)
marcetieae_coordFlagged <- marcetieae_coord %>% clean_coordinates(lon = "longitude",
                                                                  lat = "latitude",
                                                                  species = "gen_sp",
                                                                  value = "flagged",
                                                                  tests = c("equal", "gbif", 
                                                                            "institutions", 
                                                                            "outliers", "seas",
                                                                            "zeros"))

invalid_coords <- marcetieae_coord[marcetieae_coordFlagged == FALSE, ]
marcetieae_coordClean <- marcetieae_coord[marcetieae_coordFlagged  == TRUE, ]

#Binding invalid_coords to marcetieae_noCoord (5,586)
marcetieae_noCoord <- rbind(marcetieae_noCoord, invalid_coords)

rm(invalid_coords, marcetieae_coordFlagged, marcetieae_coord)

#=====================================================================================================#

#===============================#
# CLEANING MUNICIPALITIES NAMES #
#===============================#

#Dividing dataset in registers with and without values for stateprovince (859 and 4,727)
marcetieae_noState <- marcetieae_noCoord %>% filter(is.na(stateprovince))
marcetieae_state <- marcetieae_noCoord %>% filter(!is.na(stateprovince))

#NA values for municipalities (1,937)
marcetieae_na <- marcetieae_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(marcetieae_noState$county))
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
marcetieae_noState$municipality_gbif_std <- correct.mun(marcetieae_noState$municipality_gbif)
marcetieae_state$municipality_gbif_std <- correct.mun(marcetieae_state$municipality_gbif)
marcetieae_noState$county_std <- correct.mun(marcetieae_noState$county)
marcetieae_state$county_std <- correct.mun(marcetieae_state$county)

#Selecting dtb municipalities in which the campos rupestres occur
dtb_cr <- dtb[dtb$Nome_Município %in% list_mun_std, ]

#Manually checking and correcting by crossing the list assigned as 'check' with the dtb_cr dataset. 
check <- data.frame(plyr::count(marcetieae_noState$municipality_gbif_std)) 
marcetieae_noState$municipality_gbif_std[marcetieae_noState$municipality_gbif_std == 
                                           "Anaje"] <- "Anage"
marcetieae_noState$municipality_gbif_std[marcetieae_noState$municipality_gbif_std == 
                                           "Gouvea"] <- "Gouveia"
marcetieae_noState$municipality_gbif_std[marcetieae_noState$municipality_gbif_std == 
                                           "Sao Joao D'Alianca"] <- "Sao Joao D'alianca"
marcetieae_noState$municipality_gbif_std[marcetieae_noState$municipality_gbif_std 
                                         %in% c("Brasil", "Goias",
                                                "Mun?")] <- NA

check <- data.frame(plyr::count(marcetieae_noState$county_std)) 
marcetieae_noState$county_std[marcetieae_noState$county_std 
                              %in% c("Brasil", "5Km Sul da Cidade")] <- NA
marcetieae_na <- rbind(marcetieae_na, 
                       marcetieae_noState[which(is.na(marcetieae_noState$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                                  is.na(marcetieae_noState$county_std)), 1:18]) 
marcetieae_noState <- marcetieae_noState[-which(is.na(marcetieae_noState$municipality_gbif_std) & #removing registers with NA values for these two attributes from marcetieae_noState
                                                  is.na(marcetieae_noState$county_std)), ]

check <- data.frame(plyr::count(marcetieae_state$municipality_gbif_std))
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Alto Garca")] <- "Alto Garcas"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Alto Paraiso","Alto Paraiso Goias",
                                           "N de Alto Paraiso")] <- "Alto Paraiso de Goias"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Abaira   Piata"] <- "Abaira"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Anaje"] <- "Anage"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Barra de Mendes"] <- "Barra do Mendes"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Betania/Floresta"] <- "Betania"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Braslilia")] <- "Brasilia"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Caitite"] <- "Caetite"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Catugi"] <- "Catuji"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Cajazeiras Sao Jose das Piranhas"] <- "Cajazeiras"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Cidade Ecletica"] <- "Santo Antonio do Descoberto"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Conceicao do Mato de Dentro",
                                           "Conc do Mato Dentro")] <- "Conceicao do Mato Dentro"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Conselheiro da Mata")] <- "Conselheiro Mata"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Contenda do Sincora")] <- "Contendas do Sincora"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Corumba  de Goias"] <- "Corumba de Goias"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Delfinopolis (?)","Delfinopoilis")] <- "Delfinopolis"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Don Basilio")] <- "Dom Basilio"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Felixandia"] <- "Felixlandia"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Gouvea"] <- "Gouveia"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Jaboticatuba")] <- "Jaboticatubas"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Olhos D#?#Agua"] <- "Olhos D'agua"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Planaltina de Goias"] <- "Planaltina"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std 
                                       %in% c("Presidente Kubitchek")] <- "Presidente Kubitschek"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std 
                                       %in% c("Pres Soares")] <- "Presidente Soares"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Rio das Contas"] <- "Rio de Contas"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Santana do Pirapama",
                                         "Santana do Pirapama")] <- "Santana de Pirapama"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std %in% 
                                         c("Santa Rita do Jacutinga")] <- "Santa Rita de Jacutinga"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std 
                                       %in% c("Santa Terezina",
                                              "Santa Teresinha")] <- "Santa Terezinha"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Sao Tome das Letras"] <- "Sao Thome das Letras"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Sao Joao da Alianca",
                                       "Sao Joao Dalianca"] <- "Sao Joao D'alianca"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std 
                                       %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                              "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std == 
                                         "Terezina de Goias", 
                                       "Mteresina de Goias"] <- "Teresina de Goias"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std 
                                       %in% c("Varzea de Palma")] <- "Varzea da Palma"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std 
                                       %in% c("Vila Bela Santissima Trindade", "Vila Bela da Sma Trindade",
                                              "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std 
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


check <- data.frame(plyr::count(marcetieae_state$county_std))
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Alto Paraiso", "Alto Paraa-so de Goias",
                                "N de Alto Paraiso")] <- "Alto Paraiso de Goias"
marcetieae_state$county_std[marcetieae_state$county_std == 
                              "Anaje"] <- "Anage"
marcetieae_state$county_std[marcetieae_state$county_std == 
                              "Barra de Mendes"] <- "Barra do Mendes"
marcetieae_state$county_std[marcetieae_state$county_std == 
                              "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Conceicao do Mato de Dentro",
                                "Conc do Mato Dentro","Conceicao Doato Dentro")] <- "Conceicao do Mato Dentro"
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Conselheiro Matta")] <- "Conselheiro Mata"
marcetieae_state$county_std[marcetieae_state$county_std == 
                              "Cristalina Mun"] <- "Cristalina"
marcetieae_state$county_std[marcetieae_state$county_std == 
                              "Don Basilio"] <- "Dom Basilio"
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Golvea")] <- "Golveia"
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Grao Mogol Mun")] <- "Grao Mogol"
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Jaboticatuba")] <- "Jaboticatubas"
marcetieae_state$county_std[marcetieae_state$county_std == 
                              "Joaquim Fela-cio"] <- "Joaquim Felicio"
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Morro do Chapeu Mun")] <- "Morro do Chapeu"
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Planaltina de Goias")] <- "Planaltina"
marcetieae_state$county_std[marcetieae_state$county_std == 
                              "Rio das Contas",
                            "Rio de Contas Mun"] <- "Rio de Contas"
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Sao Joao da Alianca", 
                                "Sao Joao D##Alianca")] <- "Sao Joao D'alianca"
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Sao Tome das Letras ")] <- "Sao Thome das Letras "
marcetieae_state$county_std[marcetieae_state$county_std %in% 
                              c("Vila Bela de Santa-ssima Trindade", 
                                "Vila Bela de Santissima Trindade",
                                "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
marcetieae_state$county_std[marcetieae_state$county_std == 
                              "Santana de Pirapama"] <- "Santana do Pirapama"
marcetieae_state$municipality_gbif_std[marcetieae_state$municipality_gbif_std 
                                       %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                              "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
marcetieae_state$county_std[marcetieae_state$county_std == 
                              "Brasa-lia"] <- "Brasilia"
marcetieae_state$county_std[marcetieae_state$county_std 
                            %in% c("", "Chapada dos Guimaraes",
                                   "Chapada dos Veadeiros", "Chapada Gaucha",
                                   "Goias", "Minas Gerais", "No Disponible",
                                   "Pe", "Ufba Ondina", "Veadeiros", 
                                   "Minas Gerais", "Mato Grosso", "Serra de Ibitipoca",
                                   "Serra do Cabral", "Serra do Cipo", "Serra do Espinhaco",
                                   "Entre Serro E Lagoa Santa","Lavras Sao Joao Del Rey",
                                   "Chapada dos Veados")] <- NA

#Binding registers with NA values for these two attributes (4,119)
marcetieae_na <- rbind(marcetieae_na, 
                       marcetieae_state[which(is.na(marcetieae_state$municipality_gbif_std) & 
                                                is.na(marcetieae_state$county_std)), 1:18])

#Removing registers with NA values for these two attributes from marcetieae_state (3,344)
marcetieae_state <- marcetieae_state[-which(is.na(marcetieae_state$municipality_gbif_std) & 
                                             is.na(marcetieae_state$county_std)), ]

#Filtering marcetieae_noState and marcetieae_state with dtb_cr (24 and 2,062)
marcetieae_noState_filt <- marcetieae_noState %>% filter (municipality_gbif_std %in% dtb_cr$Nome_Município |
                                                            county_std %in% dtb_cr$Nome_Município)
marcetieae_state_filt <- marcetieae_state %>% filter (municipality_gbif_std %in% dtb_cr$Nome_Município |
                                                        county_std %in% dtb_cr$Nome_Município)

#=====================================================================================================#

#================================#
# INFERRING MUNICIPALITIES NAMES # 
#================================#

#Standardizing 'locality'
marcetieae_na$locality_std <- correct.mun(marcetieae_na$locality)

#Vector with municipalities names from dtb to be used with 'grepl'. This is important to match the whole municipality's name
grepl_munc_concla_cr <- c()
for(i in 1:nrow(dtb_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", dtb_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(marcetieae_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(marcetieae_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = marcetieae_na$locality_std[j])){
      vec[j] <- dtb_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(marcetieae_na$locality_std)] <- NA

#New column with inferred municipality names
marcetieae_na$municipality <- vec

#Removing registers with NA for 'municipality' (896)
marcetieae_na <- marcetieae_na %>% filter(!is.na(municipality))

#Concatenating information on municipality into a unique column 
marcetieae_noState_filt$municipality <- NA
for(i in 1:nrow(marcetieae_noState)){
  if(!is.na(marcetieae_noState_filt$municipality_gbif_std[i])){
    marcetieae_noState_filt$municipality[i] <- marcetieae_noState_filt$municipality_gbif_std[i] 
  } else if(!is.na(marcetieae_noState_filt$county_std[i])){
    marcetieae_noState_filt$municipality[i] <- marcetieae_noState_filt$county_std[i]
  }
}

marcetieae_state_filt$municipality <- NA
for(i in 1:nrow(marcetieae_state)){
  if(!is.na(marcetieae_state_filt$municipality_gbif_std[i])){
    marcetieae_state_filt$municipality[i] <- marcetieae_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(marcetieae_state_filt$county_std[i])){
    marcetieae_state_filt$municipality[i] <- marcetieae_state_filt$county_std[i]
  }
}

#Concatenating data sets without coordinates (2,982)
marcetieae_noState_filt <- marcetieae_noState_filt %>% 
  dplyr::select(colnames(marcetieae_noState_filt)[!colnames(marcetieae_noState_filt) %in%
                                                    c("county","county_std",
                                                      "municipality_gbif",
                                                      "municipality_gbif_std")])

marcetieae_state_filt <- marcetieae_state_filt %>% 
  dplyr::select(colnames(marcetieae_state_filt)[!colnames(marcetieae_state_filt) %in%
                                                  c("county","county_std",
                                                    "municipality_gbif",
                                                    "municipality_gbif_std")])

marcetieae_na <- marcetieae_na %>% 
  dplyr::select(colnames(marcetieae_na)[!colnames(marcetieae_na) %in%
                                          c("county","county_std",
                                            "municipality_gbif",
                                            "municipality_gbif_std")])

marcetieae_noCoord_inf <- rbind(marcetieae_noState_filt, marcetieae_state_filt, marcetieae_na, fill = TRUE)

#Registers occurring in homonyms municipalities (39)
reg_hom <- marcetieae_noCoord_inf[marcetieae_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
marcetieae_noCoord_inf <- marcetieae_noCoord_inf %>% filter(!municipality %in% homonyms)
marcetieae_noCoord_inf <- rbind(marcetieae_noCoord_inf, reg_hom) #binding them after removing those without state/province information

rm(marcetieae_na, marcetieae_noState, marcetieae_noState_filt, marcetieae_state, marcetieae_state_filt,
   dtb, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, marcetieae_noCoord, reg_hom, check)

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
for(i in 1:length(marcetieae_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(marcetieae_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      marcetieae_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      marcetieae_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Inferring coordinates for homonyms municipalities. Here, I relied on the combination between state and municipality. If an specific record lacks information on any of those attributes, coordinates are not inferred.
for(i in 1:length(marcetieae_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(marcetieae_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       marcetieae_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(marcetieae_noCoord_inf$stateprovince[i])){
      marcetieae_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      marcetieae_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing records without coordinates (2,952)
marcetieae_noCoord_inf <- marcetieae_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, dtb_cr, centr_hom, homonyms, not_included)

#=====================================================================================================#

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unimportant columns (7,261)
marcetieae_noCoord_inf <- marcetieae_noCoord_inf %>% dplyr::select(id,
                                                                   institutioncode,
                                                                   collectioncode,
                                                                   catalognumber,
                                                                   gen_sp,
                                                                   subspecies,
                                                                   identifiedby,
                                                                   latitude,
                                                                   longitude)

marcetieae_coordClean <- marcetieae_coordClean %>% dplyr::select(id,
                                                                 institutioncode,
                                                                 collectioncode,
                                                                 catalognumber,
                                                                 gen_sp,
                                                                 subspecies,
                                                                 identifiedby,
                                                                 latitude,
                                                                 longitude)

marcetieae_clean <- rbind(marcetieae_coordClean, marcetieae_noCoord_inf)

#Is there any duplicated registers? If so, removing them (6,848)
marcetieae_clean[duplicated(marcetieae_clean$id) == TRUE, ]
marcetieae_clean <- unique(marcetieae_clean) 

rm(marcetieae, marcetieae_coordClean, marcetieae_noCoord_inf, i, j, correct.mun, generate.names,
   replace.names, titling)

#Writing *.csv 
write.csv(marcetieae_clean, file = "datasets/Marcetieae/marcetieae_clean.csv", row.names = FALSE)

#====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Standardizing gen_sp column
marcetieae_clean$gen_sp <- gsub(" ", "_", marcetieae_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the campos rupestres shapefile and, then, with the grids (3,150)
coords <- marcetieae_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- marcetieae_clean
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
write.csv(coords, "datasets/Marcetieae/marcetieae_cr.csv", row.names = FALSE)
