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

#Reading GBIF data (16,160)
calli_gbif <- fread(file = "datasets/Calliandra/0155757-200613084148143/occurrence.txt",
                      na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Reducing data dimensionality by selecting only necessary columns
calli_gbif <- calli_gbif %>% dplyr::select(institutionCode,
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
calli_gbif <- calli_gbif %>% rename("species" = specificEpithet,
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
calli_gbif <- cbind(id = 1:nrow(calli_gbif), calli_gbif)

#=============#
# speciesLink #
#=============#

#Reading spLink (14,618)
calli_spLink <- fread(file = "datasets/Calliandra/speciesLink_all_96539_20210114211943.txt", 
                        na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Selecting important attributes
calli_spLink <- calli_spLink %>% dplyr::select(institutioncode,
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
calli_spLink$longitude <- as.numeric(as.character(calli_spLink$longitude))
calli_spLink$latitude <- as.numeric(as.character(calli_spLink$latitude))

#Giving an unique ID number for each record
calli_spLink <- cbind(id = (nrow(calli_gbif) + 1):(nrow(calli_gbif) + nrow(calli_spLink)), calli_spLink)

#======================================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (30,778), and adding a column to define the original dataset
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

calli <- merge.with.source(x = calli_gbif,
                             y = calli_spLink,
                             name.x = "gbif",
                             name.y = "splink")

rm(merge.with.source, calli_gbif, calli_spLink)

#======================================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (observations with NA for any of those attributes
#were not considered) (23,540)
a <- calli
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
calli <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
calli <- calli %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
calli$latitude[calli$latitude == 0] <- NA
calli$longitude[calli$longitude == 0] <- NA

#Removing registers without determiner's name (13,743) 
#plyr::count(calli$identifiedby)
calli$identifiedby[calli$identifiedby %in% c("?", "-", "#", "##", "//05/1987",
                                                 "24/X/2013", "`") | calli$identifiedby == 0] <- NA
calli <- calli %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- calli$stateprovince
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
calli$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
calli$stateprovince[calli$stateprovince == "?" | calli$stateprovince == "-"] <- NA
plyr::count(calli$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur (9,038)
#according to Silveira et al (2016) 
calli <- calli %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                          "Minas Gerais", 
                                                                          "Bahia", 
                                                                          "Pernambuco", 
                                                                          "Paraiba", 
                                                                          "Mato Grosso",
                                                                          "Distrito Federal"))

#Removing records without species level identification (8,772)
calli <- calli %>% filter(!is.na(species))

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (8,723)
#plyr::count(calli$basisofrecord)
calli <- calli %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                                   "s", "S", "PreservedSpecimen", "PreserverdSpecimen"))
calli <- calli %>% dplyr::select(-basisofrecord)

#=====================================================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

#Creating a vector to work with
identifiedby <- calli$identifiedby

#Extracting a vector including determiners' names. It is preferable to use this vector instead of the data set itself because any changes are easily reversible. 
identifiedby <- as.character(calli$identifiedby)

#Ideally, it is better to start with a list of taxonomic specialists' compiled beforehand. Alternatively, as experts likely indentified the majority of samples from a given taxon, it is possible to infer specialists based on identification frequency. In this example we looked for specialists at the top of the list below. 
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#To improve accuracy, we confirmed if names in 'names.count' with at least 3 identifications were specialists by searching for taxonomic publications for the family of the focal group and authored by each name at Google Scholar. In this example, we searched for: allintitle: Leguminosae OR Calliandra author:"determiner".

#Next, based on the function 'replace.names', we standardized specialist's name. This is done in two iterations:
#(1) The first iteration returns, for manual evaluation, the automatically replaced names (names above the 'top' threshold) and names that are worth to be checked (names above the 'bottom' threshold but below the 'top' threshold).
#(2) In the second iteration, names that were erroneously replaced in the first iteration should be included in the argument 'not replace'. Likewise, names that were supposed to be replaced but were below the 'top' threshold should be included in the argument 'replace'.

#Because the procedure is iterative, the user should not change the vector 'identifiedBy' directly in the first iteration. For this purpose, we used a secondary vector ('identifiedBy_2'). After the second iteration, when everything should be set, the user must assign 'identifiedBy' as 'identifiedBy_2' before running the protocol for the following name. 

#At the end of the standardizing procedure, the following vector should contain all specialists' names
specialists <- c()

#RC Barneby
replace.by <- "RC Barneby"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("R. C. Barneby & J. Grimes", "R. Barneby & J. Grimes",
                                            "in dupl.: Barneby, R. (HUEFS)","Barneby exsiccatae, ex. num. cit.",
                                            "Grimes, J. and Barneby, R. 1984", "R. Bameby", "R. C. Barneby; J. W. Grimes",
                                            "Rupert C. Barneby","R.C. Barneby & J. Grimes","Barneby, R. and Grimes, J. 1983",
                                            "R. Barneby & J. Grimes (NY) 1989"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#ER de Souza
replace.by <- "ER de Souza"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Souza, E.R. de", "Souza, ER",
                                            "Souza, E.R.", "E. R. de Sousa",
                                            "Soura, E.R.", "Souza, EB de",
                                            "E. R. de Souza (HUEFS) 2006-03",
                                            "Souza, E.","Souza","Sousa, E.R. de",
                                            "Sousa, ER de","E.R.Souza","Souza, ER de",
                                            "E. R. de Souza & L. P. de Queiroz",
                                            "Sousa, ER"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#GP Lewis
replace.by <- "GP Lewis"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("E.P.Lewis", "G. P. Lewis; B. B. Klitgaard",
                                            "E.P. LEWIS", "G. P. Lewis (K)"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#SA Harris
replace.by <- "SA Harris"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LP de Queiroz
replace.by <- "LP de Queiroz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("L. P. de Queiroz (HUNEB) 2002-01",
                                            "L.P.de Queiroz & R.M.Santos",
                                            "Queiroz, L", "Queiróz, LP", "Queiroz, L.P. DE",
                                            "Queiroz, LP","Queiroz, LP de; Lima, MPM de",
                                            "L.P.Queiroz","L.P.Queiroz; M.P.P. de Lima",
                                            "L. P. Queiroz","Queiroz, L.P. de",
                                            "Queiróz, L. P. de","L.P. Queiroz. L.P.",
                                            "Queiroz, LP de","Queiróz, L.P. de","Queiroz, L.P."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#ML Guedes
replace.by <-"ML Guedes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JM Fernandes
replace.by <- "JM Fernandes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Fernandes, A.; Bezerra, P.",
                                                "Fernandes, A.; Nunes, E.",
                                                "Fernandes, A & Nunes, E",
                                                "Fernandes, A.","Fernandes, CR",
                                                "A. Fernandes","Fernandes, M.H.",
                                                "Fernandes"),
                                replace = c("José Martins Fernandes","F.M.Fernandes"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#SA Renvoize
replace.by <- "SA Renvoize"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LM Borges
replace.by <- "LM Borges"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Leonardo M. Borges"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#S Leython
replace.by <- "S Leython"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("S Leython"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#RP Clark
replace.by <- "RP Clark"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("RP Clark"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#JW Grimes
replace.by <- "JW Grimes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("James Grimes","\"\"Grimes, J.W.; Barneby, R.C.\"\"",
                                            "J. W. Grimes; R. C. Barneby","J. W. Grimes & R. C. Barneby"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#HS Irwin
replace.by <- "HS Irwin"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("H. S. Irwin; J. W. Grimes"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#VC Souza
replace.by <- "VC de Souza"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Sousa, VF"),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#CEB Proença
replace.by <- "CEB Proença"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("C. E. B. Proença", "Proença, C.; Pereira, F."))
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#GCL Vasconcelos
replace.by <- "GCL Vasconcelos"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#A Fernandes
replace.by <- "A Fernandes"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("JM Fernandes","Fernandes, CR","Fernandes, M.H.",
                                                "Fernandes","A.Fernandes; Nunes"),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#W Mantovani
replace.by <- "W Mantovani"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#CW Fagg
replace.by <- "CW Fagg"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#HC de Lima
replace.by <- "HC de Lima"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("M. P. M. de Lima e H. C. de Lima",
                                            "Lima, H.C.","Haroldo C. de Lima",
                                            "H.C. Lima & A. Novaes","Lima, HC",
                                            "Lima, HC"))
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#A Ducke
replace.by <- "A Ducke"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LSB Jordão
replace.by <- "LSB Jordão"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#L Rico
replace.by <- "L Rico"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#MF Simon
replace.by <- "MF Simon"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#PHW Taubert
replace.by <- "PHW Taubert"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
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

#JS Silva
replace.by <- "JS Silva"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.8, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c("Silva, SDP","Silva, C.A.S. da","Silva, M.G. da",
                                                "Silva, MA", "F.F.S. Silva", "S.D.P. Silva", "Silva, DR",
                                                "Silva, J.L.", "Silva, FC", "Silva, R.A.", "Silva, F.F.S.",
                                                "Silva, S.D.P.", "Silva, N.T.", "Silva", "Silva M.A.",
                                                "Silveira A.M.", "SILVEIRA, J.E.", "Silva, M.R.",
                                                "Silva, MG da"),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#A Bocage
replace.by <- "A Bocage"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.8, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

# E Forero
replace.by <- "E Forero"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.8, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2 
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#Replacing column
calli$identifiedby <- identifiedby

#Filtering (6,779)
calli <- calli %>% filter(identifiedby %in% specialists)

rm(identifiedby, identifiedby_2, specialists, names.count, replace.by)

#======================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

library(flora)

#Generating a column for scientific names (without authors and including infraspecific epithet)
calli$gen_sp <- paste(calli$genus,
                        calli$species,
                        calli$subspecies,
                        sep = " ")

#Names
taxa <- plyr::count(calli$gen_sp)
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
#write.csv(taxa_suggested, file = "lists/Calliandra/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Calliandra/taxa_corrected.csv", stringsAsFactors = F, 
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

#Removing invalid taxa (6,564)
calli$gen_sp <- gsub("NA", "", calli$gen_sp)
calli$gen_sp <- trimws(calli$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen) | taxa_gensp$gen != "Calliandra"]
calli <- calli[!calli$gen_sp %in% invalid_taxa, ]

#Correcting the data set (varieties suppressed)
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(calli)){
  for(j in 1:nrow(taxa_gensp)){
    if(calli$gen_sp[i] == taxa_gensp$replace[j]){
      calli$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

rm(taxa, taxa_suggested, taxa_corrected, taxa_gensp, invalid_taxa, str)

#======================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (3,894)
calli_coord <- calli %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (2,670)
calli_noCoord <- calli %>% filter(is.na(latitude) | is.na(longitude))

#Removing invalid coordinates and transferring them to calli_noCoord (invalid coordinates are indicated by function clean_coordinates() in the next session)
calli_noCoord <- rbind(calli_noCoord, calli_coord[c(3152), ])
calli_coord <- calli_coord[-c(3152), ]

#======================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (3,796)
calli_coordFlagged <- calli_coord %>% clean_coordinates(lon = "longitude",
                                                            lat = "latitude",
                                                            species = "gen_sp",
                                                            value = "flagged",
                                                            tests = c("equal", "gbif", 
                                                                      "institutions", 
                                                                      "outliers", "seas",
                                                                      "zeros"))

invalid_coords <- calli_coord[calli_coordFlagged == FALSE, ]
calli_coordClean <- calli_coord[calli_coordFlagged  == TRUE, ]

#Binding invalid_coords to calli_noCoord (2,768)
calli_noCoord <- rbind(calli_noCoord, invalid_coords)

rm(invalid_coords, calli_coordFlagged, calli_coord)

#======================================================================================#

#===============================#
# CLEANING MUNICIPALITIES NAMES #
#===============================#

#Dividing dataset in registers with and without values for stateprovince (341 and 2,427)
calli_noState <- calli_noCoord %>% filter(is.na(stateprovince))
calli_state <- calli_noCoord %>% filter(!is.na(stateprovince))

#NA values for municipalities (797)
calli_na <- calli_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(calli_noState$county))
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
calli_noState$municipality_gbif_std <- correct.mun(calli_noState$municipality_gbif)
calli_state$municipality_gbif_std <- correct.mun(calli_state$municipality_gbif)
calli_noState$county_std <- correct.mun(calli_noState$county)
calli_state$county_std <- correct.mun(calli_state$county)

#Selecting dtb municipalities in which the campos rupestres occur
dtb_cr <- dtb[dtb$Nome_Município %in% list_mun_std, ]

#Manually checking and correcting by crossing the list assigned as 'check' with the dtb_cr dataset. 
#check <- data.frame(plyr::count(calli_noState$municipality_gbif_std)) #ok
calli_noState$municipality_gbif_std[calli_noState$municipality_gbif_std == 
                                        "Anaje"] <- "Anage"
calli_noState$municipality_gbif_std[calli_noState$municipality_gbif_std == 
                                        "Gouvea"] <- "Gouveia"
calli_noState$municipality_gbif_std[calli_noState$municipality_gbif_std == 
                                        "Sao Joao D'Alianca"] <- "Sao Joao D'alianca"
calli_noState$municipality_gbif_std[calli_noState$municipality_gbif_std 
                                      %in% c("Brasil", "Goias",
                                             "Mun?")] <- NA

#check <- data.frame(plyr::count(calli_noState$county_std)) #ok
calli_noState$county_std[calli_noState$county_std 
                           %in% c("Brasil", "5Km Sul da Cidade")] <- NA
calli_na <- rbind(calli_na, 
                    calli_noState[which(is.na(calli_noState$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                            is.na(calli_noState$county_std)), 1:18]) 
calli_noState <- calli_noState[-which(is.na(calli_noState$municipality_gbif_std) & #removing registers with NA values for these two attributes from calli_noState
                                            is.na(calli_noState$county_std)), ]

#check <- data.frame(plyr::count(calli_state$municipality_gbif_std))
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std %in% 
                                      c("Alto Garca")] <- "Alto Garcas"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std %in% 
                                      c("Alto Paraiso","Alto Paraiso Goias")] <- "Alto Paraiso de Goias"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Abaira   Piata"] <- "Abaira"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Anaje"] <- "Anage"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Barra de Mendes"] <- "Barra do Mendes"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Betania/Floresta"] <- "Betania"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std %in% 
                                      c("Braslilia")] <- "Brasilia"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                    "Caitite"] <- "Caetite"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Catugi"] <- "Catuji"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Cajazeiras Sao Jose das Piranhas"] <- "Cajazeiras"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Cidade Ecletica"] <- "Santo Antonio do Descoberto"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std %in% 
                                      c("Conceicao do Mato de Dentro",
                                        "Conc do Mato Dentro")] <- "Conceicao do Mato Dentro"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std %in% 
                                    c("Contenda do Sincora")] <- "Contendas do Sincora"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Corumba  de Goias"] <- "Corumba de Goias"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std %in% 
                                    c("Delfinopolis (?)","Delfinopoilis")] <- "Delfinopolis"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std %in% 
                                      c("Don Basilio")] <- "Dom Basilio"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Felixandia"] <- "Felixlandia"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Gouvea"] <- "Gouveia"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std %in% 
                                      c("Jaboticatuba")] <- "Jaboticatubas"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Olhos D#?#Agua"] <- "Olhos D'agua"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Planaltina de Goias"] <- "Planaltina"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std 
                                    %in% c("Presidente Kubitchek")] <- "Presidente Kubitschek"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Rio das Contas"] <- "Rio de Contas"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Santana do Pirapama"] <- "Santana de Pirapama"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Sao Tome das Letras"] <- "Sao Thome das Letras"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Sao Joao da Alianca"] <- "Sao Joao D'alianca"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std 
                                    %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                           "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std == 
                                      "Terezina de Goias"] <- "Teresina de Goias"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std 
                                    %in% c("Varzea de Palma")] <- "Varzea da Palma"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std 
                                    %in% c("Vila Bela Santissima Trindade", "Vila Bela da Sma Trindade",
                                           "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std 
                                    %in% c("Chapada da Cotagem", "Chapada Diamantina",
                                           "Chapada do Araripe", "Chapada dos Guimaraes",
                                           "Chapada dos Veadeiros", "Paraiba", "Pernambuco",
                                           "Serra do Cipo", "Serra do Espinhaco", "Serra do Tombador",
                                           "Ufba Ondina", "Veadeiros", "Bahia", "Chapadao do Ceu", "Goias", 
                                           "Minas Gerais", "Serra do Cipo", "Serra do Espinhaco", 
                                           "Serra do Tombador", "Veadeiros", "Bahia", "Betim/Brumadinho",
                                           "Coletada A Margem da Estrada de Macambinha","Br 135 Km 404",
                                           "Lavras Sao Joao Del Rey", "Ouro Preto/ Mariana",
                                           "Mato Grosso")] <- NA


#check <- data.frame(plyr::count(calli_state$county_std))
calli_state$county_std[calli_state$county_std %in% 
                           c("Alto Paraiso", "Alto Paraa-so de Goias")] <- "Alto Paraiso de Goias"
calli_state$county_std[calli_state$county_std == 
                           "Anaje"] <- "Anage"
calli_state$county_std[calli_state$county_std == 
                           "Barra de Mendes"] <- "Barra do Mendes"
calli_state$county_std[calli_state$county_std == 
                           "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std %in% 
                                      c("Conceicao do Mato de Dentro",
                                        "Conc do Mato Dentro","Conceicao Doato Dentro")] <- "Conceicao do Mato Dentro"
calli_state$county_std[calli_state$county_std == 
                           "Cristalina Mun"] <- "Cristalina"
calli_state$county_std[calli_state$county_std == 
                         "Don Basilio"] <- "Dom Basilio"
calli_state$county_std[calli_state$county_std %in% 
                           c("Golvea")] <- "Golveia"
calli_state$county_std[calli_state$county_std %in% 
                           c("Grao Mogol Mun")] <- "Grao Mogol"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std %in% 
                                      c("Jaboticatuba")] <- "Jaboticatubas"
calli_state$county_std[calli_state$county_std == 
                           "Joaquim Fela-cio"] <- "Joaquim Felicio"
calli_state$county_std[calli_state$county_std %in% 
                           c("Morro do Chapeu Mun")] <- "Morro do Chapeu"
calli_state$county_std[calli_state$county_std %in% 
                         c("Planaltina de Goias")] <- "Planaltina"
calli_state$county_std[calli_state$county_std == 
                           "Rio das Contas",
                       "Rio de Contas Mun"] <- "Rio de Contas"
calli_state$county_std[calli_state$county_std %in% 
                           c("Sao Joao da Alianca", 
                             "Sao Joao D##Alianca")] <- "Sao Joao D'alianca"
calli_state$county_std[calli_state$county_std %in% 
                           c("Sao Tome das Letras ")] <- "Sao Thome das Letras "
calli_state$county_std[calli_state$county_std %in% 
                           c("Vila Bela de Santa-ssima Trindade", 
                             "Vila Bela de Santissima Trindade",
                             "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
calli_state$county_std[calli_state$county_std == 
                           "Santana de Pirapama"] <- "Santana do Pirapama"
calli_state$municipality_gbif_std[calli_state$municipality_gbif_std 
                                    %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                           "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
calli_state$county_std[calli_state$county_std == 
                           "Brasa-lia"] <- "Brasilia"
calli_state$county_std[calli_state$county_std 
                         %in% c("", "Chapada dos Guimaraes",
                                "Chapada dos Veadeiros", "Chapada Gaucha",
                                "Goias", "Minas Gerais", "No Disponible",
                                "Pe", "Ufba Ondina", "Veadeiros", 
                                "Minas Gerais", "Mato Grosso", "Serra de Ibitipoca",
                                "Serra do Cabral", "Serra do Cipo", "Serra do Espinhaco",
                                "Entre Serro E Lagoa Santa","Lavras Sao Joao Del Rey")] <- NA

#Binding registers with NA values for these two attributes (1,661)
calli_na <- rbind(calli_na, 
                    calli_state[which(is.na(calli_state$municipality_gbif_std) & 
                                          is.na(calli_state$county_std)), 1:18])

#Removing registers with NA values for these two attributes from calli_state (1,266)
calli_state <- calli_state[which(!is.na(calli_state$municipality_gbif_std) & 
                                        is.na(calli_state$county_std)), ]

#Filtering calli_noState and calli_state with dtb_cr (96 and 941)
calli_noState_filt <- calli_noState %>% filter (municipality_gbif_std %in% dtb_cr$Nome_Município |
                                                      county_std %in% dtb_cr$Nome_Município)
calli_state_filt <- calli_state %>% filter (municipality_gbif_std %in% dtb_cr$Nome_Município |
                                                  county_std %in% dtb_cr$Nome_Município)

#=====================================================================================================#

#================================#
# INFERRING MUNICIPALITIES NAMES # 
#================================#

#Standardizing 'locality'
calli_na$locality_std <- correct.mun(calli_na$locality)

#Vector with municipalities names from dtb to be used with 'grepl'. This is important to match the whole municipality's name
grepl_munc_concla_cr <- c()
for(i in 1:nrow(dtb_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", dtb_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(calli_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(calli_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = calli_na$locality_std[j])){
      vec[j] <- dtb_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(calli_na$locality_std)] <- NA

#New column with inferred municipality names
calli_na$municipality <- vec

#Removing registers with NA for 'municipality'
calli_na <- calli_na %>% filter(!is.na(municipality))

#Concatenating information on municipality into a unique column 
calli_noState_filt$municipality <- NA
for(i in 1:nrow(calli_noState)){
  if(!is.na(calli_noState_filt$municipality_gbif_std[i])){
    calli_noState_filt$municipality[i] <- calli_noState_filt$municipality_gbif_std[i] 
  } else if(!is.na(calli_noState_filt$county_std[i])){
    calli_noState_filt$municipality[i] <- calli_noState_filt$county_std[i]
  }
}

calli_state_filt$municipality <- NA
for(i in 1:nrow(calli_state)){
  if(!is.na(calli_state_filt$municipality_gbif_std[i])){
    calli_state_filt$municipality[i] <- calli_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(calli_state_filt$county_std[i])){
    calli_state_filt$municipality[i] <- calli_state_filt$county_std[i]
  }
}

#Concatenating data sets without coordinates (1,667)
calli_noState_filt <- calli_noState_filt %>% 
  dplyr::select(colnames(calli_noState_filt)[!colnames(calli_noState_filt) %in%
                                                 c("county","county_std",
                                                   "municipality_gbif",
                                                   "municipality_gbif_std")])

calli_state_filt <- calli_state_filt %>% 
  dplyr::select(colnames(calli_state_filt)[!colnames(calli_state_filt) %in%
                                               c("county","county_std",
                                                 "municipality_gbif",
                                                 "municipality_gbif_std")])

calli_na <- calli_na %>% 
  dplyr::select(colnames(calli_na)[!colnames(calli_na) %in%
                                       c("county","county_std",
                                         "municipality_gbif",
                                         "municipality_gbif_std")])

calli_noCoord_inf <- rbind(calli_noState_filt, calli_state_filt, calli_na, fill = TRUE)

#Registers occurring in homonyms municipalities (1,677)
reg_hom <- calli_noCoord_inf[calli_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
calli_noCoord_inf <- calli_noCoord_inf %>% filter(!municipality %in% homonyms)
calli_noCoord_inf <- rbind(calli_noCoord_inf, reg_hom) #binding them after removing those without state/province information

rm(calli_na, calli_noState, calli_noState_filt, calli_state, calli_state_filt,
   dtb, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, calli_noCoord, reg_hom)

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
for(i in 1:length(calli_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(calli_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      calli_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      calli_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Inferring coordinates for homonyms municipalities. Here, I relied on the combination between state and municipality. If an specific record lacks information on any of those attributes, coordinates are not inferred.
for(i in 1:length(calli_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(calli_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       calli_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(calli_noCoord_inf$stateprovince[i])){
      calli_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      calli_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing records without coordinates (811)
calli_noCoord_inf <- calli_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, dtb_cr, centr_hom, homonyms)

#=================================================================================================#

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unimportant columns (4,607)
calli_noCoord_inf <- calli_noCoord_inf %>% dplyr::select(id,
                                                             institutioncode,
                                                             collectioncode,
                                                             catalognumber,
                                                             gen_sp,
                                                             subspecies,
                                                             identifiedby,
                                                             latitude,
                                                             longitude)

calli_coordClean <- calli_coordClean %>% dplyr::select(id,
                                                           institutioncode,
                                                           collectioncode,
                                                           catalognumber,
                                                           gen_sp,
                                                           subspecies,
                                                           identifiedby,
                                                           latitude,
                                                           longitude)

calli_clean <- rbind(calli_coordClean, calli_noCoord_inf)

#Is there any duplicated registers? If so, removing them (4,455)
calli_clean[duplicated(calli_clean$id) == TRUE, ]
calli_clean <- unique(calli_clean) 

rm(calli, calli_coordClean, calli_noCoord_inf, i, j, correct.mun, generate.names,
   replace.names, titling, not_included)

#Writing *.csv 
write.csv(calli_clean, file = "datasets/Calliandra/calli_clean.csv", row.names = FALSE)

#====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Standardizing gen_sp column
calli_clean$gen_sp <- gsub(" ", "_", calli_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the campos rupestres shapefile and, then, with the grids (1,686)
coords <- calli_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- calli_clean
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
write.csv(coords, "datasets/Calliandra/calli_cr.csv", row.names = FALSE)
