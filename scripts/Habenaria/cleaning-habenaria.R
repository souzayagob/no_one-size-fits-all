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

#Reading GBIF data (14,677)
habenaria_gbif <- fread(file = "datasets/Habenaria/0155769-200613084148143/occurrence.txt",
                   na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Reducing data dimensionality by selecting only necessary columns
habenaria_gbif <- habenaria_gbif %>% dplyr::select(institutionCode,
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
habenaria_gbif <- habenaria_gbif %>% rename("species" = specificEpithet,
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
habenaria_gbif <- cbind(id = 1:nrow(habenaria_gbif), habenaria_gbif)

#=============#
# speciesLink #
#=============#

#Reading spLink (12,603)
habenaria_spLink <- fread(file = "datasets/Habenaria/speciesLink_all_105061_20210114213510.txt", 
                     na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")

#Selecting important attributes
habenaria_spLink <- habenaria_spLink %>% dplyr::select(institutioncode,
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
habenaria_spLink$longitude <- as.numeric(as.character(habenaria_spLink$longitude))
habenaria_spLink$latitude <- as.numeric(as.character(habenaria_spLink$latitude))

#Giving an unique ID number for each record
habenaria_spLink <- cbind(id = (nrow(habenaria_gbif) + 1):(nrow(habenaria_gbif) + nrow(habenaria_spLink)), habenaria_spLink)

#======================================================================================#

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

habenaria <- merge.with.source(x = habenaria_gbif,
                          y = habenaria_spLink,
                          name.x = "gbif",
                          name.y = "splink")

rm(merge.with.source, habenaria_gbif, habenaria_spLink)

#======================================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (observations with NA for any of those attributes
#were not considered) (21,672)
a <- habenaria
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
habenaria <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
habenaria <- habenaria %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
habenaria$latitude[habenaria$latitude == 0] <- NA
habenaria$longitude[habenaria$longitude == 0] <- NA

#Removing registers without determiner's name (3,620) 
#plyr::count(habenaria$identifiedby)
habenaria$identifiedby[habenaria$identifiedby %in% c("?", "-", "#", "##", "//05/1987",
                                           "24/X/2013", "`") | habenaria$identifiedby == 0] <- NA
habenaria <- habenaria %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- habenaria$stateprovince
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
habenaria$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
habenaria$stateprovince[habenaria$stateprovince == "?" | habenaria$stateprovince == "-"] <- NA
plyr::count(habenaria$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur (7,561)
#according to Silveira et al (2016) 
habenaria <- habenaria %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                    "Minas Gerais", 
                                                                    "Bahia", 
                                                                    "Pernambuco", 
                                                                    "Paraiba", 
                                                                    "Mato Grosso",
                                                                    "Distrito Federal"))

#Removing records without species level identification (7,216)
habenaria <- habenaria %>% filter(!is.na(species))

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (7,195)
plyr::count(habenaria$basisofrecord)
habenaria <- habenaria %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                             "s", "S", "PreservedSpecimen", "PreserverdSpecimen"))
habenaria <- habenaria %>% dplyr::select(-basisofrecord)

#=====================================================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

#Creating a vector to work with
identifiedby <- habenaria$identifiedby

#Extracting a vector including determiners' names. It is preferable to use this vector instead of the data set itself because any changes are easily reversible. 
identifiedby <- as.character(habenaria$identifiedby)

#Ideally, it is better to start with a list of taxonomic specialists' compiled beforehand. Alternatively, as experts likely indentified the majority of samples from a given taxon, it is possible to infer specialists based on identification frequency. In this example we looked for specialists at the top of the list below. 
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#To improve accuracy, we confirmed if names in 'names.count' with at least 3 identifications were specialists by searching for taxonomic publications for the family of the focal group and authored by each name at Google Scholar. In this example, we searched for: allintitle: Orchidaceae OR Habenaria author:"determiner".

#Next, based on the function 'replace.names', we standardized specialist's name. This is done in two iterations:
#(1) The first iteration returns, for manual evaluation, the automatically replaced names (names above the 'top' threshold) and names that are worth to be checked (names above the 'bottom' threshold but below the 'top' threshold).
#(2) In the second iteration, names that were erroneously replaced in the first iteration should be included in the argument 'not replace'. Likewise, names that were supposed to be replaced but were below the 'top' threshold should be included in the argument 'replace'.

#Because the procedure is iterative, the user should not change the vector 'identifiedBy' directly in the first iteration. For this purpose, we used a secondary vector ('identifiedBy_2'). After the second iteration, when everything should be set, the user must assign 'identifiedBy' as 'identifiedBy_2' before running the protocol for the following name. 

#At the end of the standardizing procedure, the following vector should contain all specialists' names
specialists <- c()

#JAN Batista
replace.by <- "JAN Batista"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("A.J.N. Batista et al.", "J. A. N. Batista 2013-02-06",
                                            "J. A. N. Batista 2013-02-04","J. A. N. Batista & L. B. Bianchetti",
                                            "João Aguiar N. Batista & A. J. Ramalho", "(BHCB) Batista, JAN",
                                            "(BHCB) Batista, JAN", "J.A.N. Batista & B.M. Carvalho",
                                            "J.A.N.Batista & L.B.Bianchetti", "J.Batsita", 
                                            "J.A.N. Bautista 2013-01-09", "J.Batista & R.C.Mota", 
                                            "João Aguiar Nogueira Batista; J.A.N.Bat.", "J. A. N. Batista 2013-03-25",
                                            "J. A. N. Batista; L. de Bem Bianchetti", "J. A. N. Batista & A. J. Ramalho",
                                            "J.A.N.Batista & Bianchetti", "B.M.Carvalho & J.A.N.Batista","J.A.N. Batista & L.B. Bianchetti",
                                            "J. H. A. Batista", "J.A.N.BATISTA & L.BIANCHETTI", "J.A.N. Bautista 2013-03-25",
                                            "J.A.N. Bautista 2013-03-25","J. A. N. Batista/Dez.1998","J. A. N. Batista & Osvaldo Neto",
                                            "J.A.N. Batista & R. Belisário", "João Aguiar Nogueira Batista",
                                            "João Aguiar Nogueira Batista; Universidade Federal de Minas Gerais; Depto. Botânica #?# ICB; J.A.N"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#GFJ Pabst
replace.by <- "GFJ Pabst"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("G. Pabst, 19 Apr","Rev. Pabst, 19 Apr","C. Pabst"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LB Bianchetti
replace.by <- "LB Bianchetti"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("L. B. Bianchetti & J. A. N. Batista","L. de B. Bianchetti",
                                            "L.B. Bianchetti, C. Maury & F.T. Andrade"))
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

#F Barros
replace.by <- "F Barros"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("F.F.V.A. Barberena & F. Barros", "Barros, F.; Vinhos, F.; Rodrigues, VT (LEFB)"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#L Menini Neto
replace.by <- "L Menini Neto"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("L.M.Neto",))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#BM Carvalho
replace.by <- "BM Carvalho"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("B. M. Carvalho & J.A.N. Pastore","B. M. Carvalho &J. A. N. Batista",
                                            "B. M. Carvalho & J. A. N. Batista"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#A Cogniaux
replace.by <- "A Cogniaux"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#LS Leoni
replace.by <- "LS Leoni"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("L.S. Leoni; P.I.S. Braga","L.S. Leoni; W. Forster"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#C van den Berg
replace.by <- "C van den Berg"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Berg, C. van den","C. van den Berg & T. E. C. Meneguzzo 2014-07-23",
                                            "Berg, C. van den; Azevedo, C.O."))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#NL Abreu
replace.by <- "NL Abreu"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#RB Singer
replace.by <- "RB Singer"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("C. R. Buzatto & R. B. Singer"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#J Renz
replace.by <- "J Renz"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#TL Vieira
replace.by <- "TL Vieira"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#EA Christenson
replace.by <- "EA Christenson"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#ALV Toscano de Brito
replace.by <- "ALV Toscano de Brito"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("A. Toscano","A. T. de Brito","A.T.de Brito",
                                            "A. L. Toscano de Brito","de Brito, T.",
                                            "Toscano, Mar"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#F Kränzlin
replace.by <- "F Kränzlin"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("Kraenzlin","Kãnzlin"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#E Pessoa
replace.by <- "E Pessoa"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#AA Vale
replace.by <- "AA Vale"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#P Cribb
replace.by <- "P Cribb"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#FFVA Barberena
replace.by <- "FFVA Barberena"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c("F.F.V.A. Barberena & J.A.N. Batista"))
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#HG Reichenbach
replace.by <- "HG Reichenbach"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#J Ribeiro
replace.by <- "J Ribeiro"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#CN Fraga
replace.by <- "CN Fraga"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#AC Brade
replace.by <- "AC Brade"
identifiedby_2 <- replace.names(x = identifiedby, top = 0.85, bottom = 0.6, 
                                check.by = generate.names(str_split(replace.by, pattern = " ", n = 2)[[1]][1],
                                                          str_split(replace.by, pattern = " ", n = 2)[[1]][2]),
                                replace.by = replace.by,
                                not.replace = c(),
                                replace = c())
identifiedby <- identifiedby_2
specialists <- c(specialists, replace.by)
names.count <- as.data.frame(plyr::count(identifiedby))[order(-as.data.frame(plyr::count(identifiedby))$freq), ]

#DN Carvalho
replace.by <- "DN Carvalho"
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
habenaria$identifiedby <- identifiedby

#Filtering (6,205)
habenaria <- habenaria %>% filter(identifiedby %in% specialists)

rm(identifiedby, identifiedby_2, specialists, names.count, replace.by)

#======================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

library(flora)

#Generating a column for scientific names (without authors and including infraspecific epithet)
habenaria$gen_sp <- paste(habenaria$genus,
                     habenaria$species,
                     habenaria$subspecies,
                     sep = " ")

#Names
taxa <- plyr::count(habenaria$gen_sp)
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
write.csv(taxa_suggested, file = "lists/Habenaria/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Habenaria/taxa_corrected.csv", stringsAsFactors = F, 
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

#Removing invalid taxa (6,178)
habenaria$gen_sp <- gsub("NA", "", habenaria$gen_sp)
habenaria$gen_sp <- trimws(habenaria$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen) | taxa_gensp$gen != "Habenaria"]
habenaria <- habenaria[!habenaria$gen_sp %in% invalid_taxa, ]

#Correcting the data set (varieties suppressed)
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(habenaria)){
  for(j in 1:nrow(taxa_gensp)){
    if(habenaria$gen_sp[i] == taxa_gensp$replace[j]){
      habenaria$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

rm(taxa, taxa_suggested, taxa_corrected, taxa_gensp, invalid_taxa, str)

#======================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (1,746)
habenaria_coord <- habenaria %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (4,432)
habenaria_noCoord <- habenaria %>% filter(is.na(latitude) | is.na(longitude))

#Removing invalid coordinates and transferring them to habenaria_noCoord (invalid coordinates are indicated by function clean_coordinates() in the next session)
habenaria_noCoord <- rbind(habenaria_noCoord, habenaria_coord[c(), ])
habenaria_coord <- habenaria_coord[-c(), ]

#======================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (1,667)
habenaria_coordFlagged <- habenaria_coord %>% clean_coordinates(lon = "longitude",
                                                      lat = "latitude",
                                                      species = "gen_sp",
                                                      value = "flagged",
                                                      tests = c("equal", "gbif", 
                                                                "institutions", 
                                                                "outliers", "seas",
                                                                "zeros"))

invalid_coords <- habenaria_coord[habenaria_coordFlagged == FALSE, ]
habenaria_coordClean <- habenaria_coord[habenaria_coordFlagged  == TRUE, ]

#Binding invalid_coords to habenaria_noCoord (4,511)
habenaria_noCoord <- rbind(habenaria_noCoord, invalid_coords)

rm(invalid_coords, habenaria_coordFlagged, habenaria_coord)

#======================================================================================#

#===============================#
# CLEANING MUNICIPALITIES NAMES #
#===============================#

#Dividing dataset in registers with and without values for stateprovince (922 and 3,589)
habenaria_noState <- habenaria_noCoord %>% filter(is.na(stateprovince))
habenaria_state <- habenaria_noCoord %>% filter(!is.na(stateprovince))

#NA values for municipalities (1,335)
habenaria_na <- habenaria_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(habenaria_noState$county))
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
habenaria_noState$municipality_gbif_std <- correct.mun(habenaria_noState$municipality_gbif)
habenaria_state$municipality_gbif_std <- correct.mun(habenaria_state$municipality_gbif)
habenaria_noState$county_std <- correct.mun(habenaria_noState$county)
habenaria_state$county_std <- correct.mun(habenaria_state$county)

#Selecting dtb municipalities in which the campos rupestres occur
dtb_cr <- dtb[dtb$Nome_Município %in% list_mun_std, ]

#Manually checking and correcting by crossing the list assigned as 'check' with the dtb_cr dataset. 
check <- data.frame(plyr::count(habenaria_noState$municipality_gbif_std)) 
habenaria_noState$municipality_gbif_std[habenaria_noState$municipality_gbif_std == 
                                     "Anaje"] <- "Anage"
habenaria_noState$municipality_gbif_std[habenaria_noState$municipality_gbif_std == 
                                     "Gouvea"] <- "Gouveia"
habenaria_noState$municipality_gbif_std[habenaria_noState$municipality_gbif_std == 
                                     "Sao Joao D'Alianca"] <- "Sao Joao D'alianca"
habenaria_noState$municipality_gbif_std[habenaria_noState$municipality_gbif_std 
                                   %in% c("Brasil", "Goias",
                                          "Mun?")] <- NA

check <- data.frame(plyr::count(habenaria_noState$county_std)) 
habenaria_noState$county_std[habenaria_noState$county_std 
                        %in% c("Brasil", "5Km Sul da Cidade")] <- NA
habenaria_na <- rbind(habenaria_na, 
                 habenaria_noState[which(is.na(habenaria_noState$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                      is.na(habenaria_noState$county_std)), 1:18]) 
habenaria_noState <- habenaria_noState[-which(is.na(habenaria_noState$municipality_gbif_std) & #removing registers with NA values for these two attributes from habenaria_noState
                                      is.na(habenaria_noState$county_std)), ]

check <- data.frame(plyr::count(habenaria_state$municipality_gbif_std))
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std %in% 
                                   c("Alto Garca")] <- "Alto Garcas"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std %in% 
                                   c("Alto Paraiso","Alto Paraiso Goias",
                                     "N de Alto Paraiso")] <- "Alto Paraiso de Goias"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Abaira   Piata"] <- "Abaira"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Anaje"] <- "Anage"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Barra de Mendes"] <- "Barra do Mendes"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Betania/Floresta"] <- "Betania"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std %in% 
                                   c("Braslilia")] <- "Brasilia"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Caitite"] <- "Caetite"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Catugi"] <- "Catuji"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Cajazeiras Sao Jose das Piranhas"] <- "Cajazeiras"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Cidade Ecletica"] <- "Santo Antonio do Descoberto"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std %in% 
                                   c("Conceicao do Mato de Dentro",
                                     "Conc do Mato Dentro")] <- "Conceicao do Mato Dentro"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std %in% 
                                   c("Conselheiro da Mata")] <- "Conselheiro Mata"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std %in% 
                                   c("Contenda do Sincora")] <- "Contendas do Sincora"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Corumba  de Goias"] <- "Corumba de Goias"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std %in% 
                                   c("Delfinopolis (?)","Delfinopoilis")] <- "Delfinopolis"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std %in% 
                                   c("Don Basilio")] <- "Dom Basilio"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Felixandia"] <- "Felixlandia"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Gouvea"] <- "Gouveia"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std %in% 
                                   c("Jaboticatuba")] <- "Jaboticatubas"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Olhos D#?#Agua"] <- "Olhos D'agua"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Planaltina de Goias"] <- "Planaltina"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std 
                                 %in% c("Presidente Kubitchek")] <- "Presidente Kubitschek"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Rio das Contas"] <- "Rio de Contas"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Santana do Pirapama"] <- "Santana de Pirapama"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Sao Tome das Letras"] <- "Sao Thome das Letras"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Sao Joao da Alianca"] <- "Sao Joao D'alianca"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std 
                                 %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                        "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std == 
                                   "Terezina de Goias", 
                                   "Mteresina de Goias"] <- "Teresina de Goias"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std 
                                 %in% c("Varzea de Palma")] <- "Varzea da Palma"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std 
                                 %in% c("Vila Bela Santissima Trindade", "Vila Bela da Sma Trindade",
                                        "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std 
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


check <- data.frame(plyr::count(habenaria_state$county_std))
habenaria_state$county_std[habenaria_state$county_std %in% 
                        c("Alto Paraiso", "Alto Paraa-so de Goias",
                          "N de Alto Paraiso")] <- "Alto Paraiso de Goias"
habenaria_state$county_std[habenaria_state$county_std == 
                        "Anaje"] <- "Anage"
habenaria_state$county_std[habenaria_state$county_std == 
                        "Barra de Mendes"] <- "Barra do Mendes"
habenaria_state$county_std[habenaria_state$county_std == 
                        "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
habenaria_state$county_std[habenaria_state$county_std %in% 
                                   c("Conceicao do Mato de Dentro",
                                     "Conc do Mato Dentro","Conceicao Doato Dentro")] <- "Conceicao do Mato Dentro"
habenaria_state$county_std[habenaria_state$county_std == 
                        "Cristalina Mun"] <- "Cristalina"
habenaria_state$county_std[habenaria_state$county_std == 
                        "Don Basilio"] <- "Dom Basilio"
habenaria_state$county_std[habenaria_state$county_std %in% 
                        c("Golvea")] <- "Golveia"
habenaria_state$county_std[habenaria_state$county_std %in% 
                        c("Grao Mogol Mun")] <- "Grao Mogol"
habenaria_state$county_std[habenaria_state$county_std %in% 
                                   c("Jaboticatuba")] <- "Jaboticatubas"
habenaria_state$county_std[habenaria_state$county_std == 
                        "Joaquim Fela-cio"] <- "Joaquim Felicio"
habenaria_state$county_std[habenaria_state$county_std %in% 
                        c("Morro do Chapeu Mun")] <- "Morro do Chapeu"
habenaria_state$county_std[habenaria_state$county_std %in% 
                        c("Planaltina de Goias")] <- "Planaltina"
habenaria_state$county_std[habenaria_state$county_std == 
                        "Rio das Contas",
                      "Rio de Contas Mun"] <- "Rio de Contas"
habenaria_state$county_std[habenaria_state$county_std %in% 
                        c("Sao Joao da Alianca", 
                          "Sao Joao D##Alianca")] <- "Sao Joao D'alianca"
habenaria_state$county_std[habenaria_state$county_std %in% 
                        c("Sao Tome das Letras ")] <- "Sao Thome das Letras "
habenaria_state$county_std[habenaria_state$county_std %in% 
                        c("Vila Bela de Santa-ssima Trindade", 
                          "Vila Bela de Santissima Trindade",
                          "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
habenaria_state$county_std[habenaria_state$county_std == 
                        "Santana de Pirapama"] <- "Santana do Pirapama"
habenaria_state$municipality_gbif_std[habenaria_state$municipality_gbif_std 
                                 %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei",
                                        "Sao Joao Del Rey")] <- "Sao Joao Del Rei"
habenaria_state$county_std[habenaria_state$county_std == 
                        "Brasa-lia"] <- "Brasilia"
habenaria_state$county_std[habenaria_state$county_std 
                      %in% c("", "Chapada dos Guimaraes",
                             "Chapada dos Veadeiros", "Chapada Gaucha",
                             "Goias", "Minas Gerais", "No Disponible",
                             "Pe", "Ufba Ondina", "Veadeiros", 
                             "Minas Gerais", "Mato Grosso", "Serra de Ibitipoca",
                             "Serra do Cabral", "Serra do Cipo", "Serra do Espinhaco",
                             "Entre Serro E Lagoa Santa","Lavras Sao Joao Del Rey",
                             "Chapada dos Veados")] <- NA

#Binding registers with NA values for these two attributes (2,791)
habenaria_na <- rbind(habenaria_na, 
                 habenaria_state[which(is.na(habenaria_state$municipality_gbif_std) & 
                                    is.na(habenaria_state$county_std)), 1:18])

#Removing registers with NA values for these two attributes from habenaria_state (2,438)
habenaria_state <- habenaria_state[which(!is.na(habenaria_state$municipality_gbif_std) & 
                                 is.na(habenaria_state$county_std)), ]

#Filtering habenaria_noState and habenaria_state with dtb_cr (29 and 1,762)
habenaria_noState_filt <- habenaria_noState %>% filter (municipality_gbif_std %in% dtb_cr$Nome_Município |
                                                county_std %in% dtb_cr$Nome_Município)
habenaria_state_filt <- habenaria_state %>% filter (municipality_gbif_std %in% dtb_cr$Nome_Município |
                                            county_std %in% dtb_cr$Nome_Município)

#=====================================================================================================#

#================================#
# INFERRING MUNICIPALITIES NAMES # 
#================================#

#Standardizing 'locality'
habenaria_na$locality_std <- correct.mun(habenaria_na$locality)

#Vector with municipalities names from dtb to be used with 'grepl'. This is important to match the whole municipality's name
grepl_munc_concla_cr <- c()
for(i in 1:nrow(dtb_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", dtb_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(habenaria_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(habenaria_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = habenaria_na$locality_std[j])){
      vec[j] <- dtb_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(habenaria_na$locality_std)] <- NA

#New column with inferred municipality names
habenaria_na$municipality <- vec

#Removing registers with NA for 'municipality' (530)
habenaria_na <- habenaria_na %>% filter(!is.na(municipality))

#Concatenating information on municipality into a unique column 
habenaria_noState_filt$municipality <- NA
for(i in 1:nrow(habenaria_noState)){
  if(!is.na(habenaria_noState_filt$municipality_gbif_std[i])){
    habenaria_noState_filt$municipality[i] <- habenaria_noState_filt$municipality_gbif_std[i] 
  } else if(!is.na(habenaria_noState_filt$county_std[i])){
    habenaria_noState_filt$municipality[i] <- habenaria_noState_filt$county_std[i]
  }
}

habenaria_state_filt$municipality <- NA
for(i in 1:nrow(habenaria_state)){
  if(!is.na(habenaria_state_filt$municipality_gbif_std[i])){
    habenaria_state_filt$municipality[i] <- habenaria_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(habenaria_state_filt$county_std[i])){
    habenaria_state_filt$municipality[i] <- habenaria_state_filt$county_std[i]
  }
}

#Concatenating data sets without coordinates (2,321)
habenaria_noState_filt <- habenaria_noState_filt %>% 
  dplyr::select(colnames(habenaria_noState_filt)[!colnames(habenaria_noState_filt) %in%
                                              c("county","county_std",
                                                "municipality_gbif",
                                                "municipality_gbif_std")])

habenaria_state_filt <- habenaria_state_filt %>% 
  dplyr::select(colnames(habenaria_state_filt)[!colnames(habenaria_state_filt) %in%
                                            c("county","county_std",
                                              "municipality_gbif",
                                              "municipality_gbif_std")])

habenaria_na <- habenaria_na %>% 
  dplyr::select(colnames(habenaria_na)[!colnames(habenaria_na) %in%
                                    c("county","county_std",
                                      "municipality_gbif",
                                      "municipality_gbif_std")])

habenaria_noCoord_inf <- rbind(habenaria_noState_filt, habenaria_state_filt, habenaria_na, fill = TRUE)

#Registers occurring in homonyms municipalities (55)
reg_hom <- habenaria_noCoord_inf[habenaria_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
habenaria_noCoord_inf <- habenaria_noCoord_inf %>% filter(!municipality %in% homonyms)
habenaria_noCoord_inf <- rbind(habenaria_noCoord_inf, reg_hom) #binding them after removing those without state/province information

rm(habenaria_na, habenaria_noState, habenaria_noState_filt, habenaria_state, habenaria_state_filt,
   dtb, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, habenaria_noCoord, reg_hom, check)

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
for(i in 1:length(habenaria_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(habenaria_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      habenaria_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      habenaria_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Inferring coordinates for homonyms municipalities. Here, I relied on the combination between state and municipality. If an specific record lacks information on any of those attributes, coordinates are not inferred.
for(i in 1:length(habenaria_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(habenaria_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       habenaria_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(habenaria_noCoord_inf$stateprovince[i])){
      habenaria_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      habenaria_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing records without coordinates (587)
habenaria_noCoord_inf <- habenaria_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, dtb_cr, centr_hom, homonyms, not_included)

#=================================================================================================#

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unimportant columns (2,254)
habenaria_noCoord_inf <- habenaria_noCoord_inf %>% dplyr::select(id,
                                                       institutioncode,
                                                       collectioncode,
                                                       catalognumber,
                                                       gen_sp,
                                                       subspecies,
                                                       identifiedby,
                                                       latitude,
                                                       longitude)

habenaria_coordClean <- habenaria_coordClean %>% dplyr::select(id,
                                                     institutioncode,
                                                     collectioncode,
                                                     catalognumber,
                                                     gen_sp,
                                                     subspecies,
                                                     identifiedby,
                                                     latitude,
                                                     longitude)

habenaria_clean <- rbind(habenaria_coordClean, habenaria_noCoord_inf)

#Is there any duplicated registers? If so, removing them (2,166)
habenaria_clean[duplicated(habenaria_clean$id) == TRUE, ]
habenaria_clean <- unique(habenaria_clean) 

rm(habenaria, habenaria_coordClean, habenaria_noCoord_inf, i, j, correct.mun, generate.names,
   replace.names, titling)

#====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Standardizing gen_sp column
habenaria_clean$gen_sp <- gsub(" ", "_", habenaria_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the campos rupestres shapefile and, then, with the grids (1,686)
coords <- habenaria_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- habenaria_clean
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
write.csv(coords, "datasets/Habenaria/habenaria_cr.csv", row.names = FALSE)