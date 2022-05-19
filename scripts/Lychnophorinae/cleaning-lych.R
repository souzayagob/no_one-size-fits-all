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
lych_gbif <- fread(file = "Datasets/Lychnophorinae/0037335-200613084148143-gbif_lychnophorinae/occurrence.txt",
                  na.strings = c("", NA), stringsAsFactors = FALSE, encoding = "UTF-8")

#Selecting important attributes
lych_gbif <- lych_gbif %>% select(institutionCode,
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
lych_gbif <- lych_gbif %>% rename("species" = specificEpithet,
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
lych_gbif <- cbind(id = 1:nrow(lych_gbif), lych_gbif)

#Encoding strings as UTF-8 (R version 4)
#class(lych_gbif) <- "data.frame" #if the current class is "data.table" "data.frame",
#the following code does not work properly. 
#for(i in 1:ncol(lych_gbif)){
#if(is.character(lych_gbif[ , i])){
#   Encoding(lych_gbif[ , i]) <- "UTF-8" 
#  }
#}

#=============#
# speciesLink #
#=============#

#Reading spLink 
lych_spLink <- fread(file = "Datasets/Lychnophorinae/speciesLink_lychnophorinae_76349_20200812120108.txt", 
                    na.strings = c("", NA), stringsAsFactors = FALSE, encoding = "UTF-8")

#Selecting important attributes
lych_spLink <- lych_spLink %>% select(institutioncode,
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
lych_spLink$longitude <- as.numeric(as.character(lych_spLink$longitude))
lych_spLink$latitude <- as.numeric(as.character(lych_spLink$latitude))

#Giving an unique ID number for each record
lych_spLink <- cbind(id = (nrow(lych_gbif) + 1):(nrow(lych_gbif) + nrow(lych_spLink)), lych_spLink)

#Encoding strings as UTF-8 
#class(lych_spLink) <- "data.frame" #if the current class is "data.table" "data.frame",
#the following code does not work properly. 
#for(i in 1:ncol(lych_spLink)){
#if(is.character(lych_spLink[ , i])){
# Encoding(lych_spLink[ , i]) <- "UTF-8" 
#}
#}

#==============================================================================#

#========================#
# CONCATENATING DATASETS #
#========================#

#Merging gbif and splink (53,613), and adding a column to define the original dataset
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

lych <- merge.with.source(x = lych_gbif,
                         y = lych_spLink,
                         name.x = "gbif",
                         name.y = "splink")

rm(merge.with.source, lych_gbif, lych_spLink)

#==============================================================================#

#================#
# PRE-REFINEMENT #
#================#

#Removing duplicated registers given the institution code, the collection code and
#the catalog number (observations with NA for any of those attributes
#were not considered) (41,599)
a <- lych
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
lych <- a
rm(a, a.prime, a.na)

#Indicating dataset specific columns
lych <- lych %>% rename("municipality_gbif" = municipality)

#Replacing 0 by NA in the coordinates 
lych$latitude[lych$latitude == 0] <- NA
lych$longitude[lych$longitude == 0] <- NA

#Removing registers without identifier name (25,684)
plyr::count(lych$identifiedby)
lych$identifiedby[lych$identifiedby %in% c("?", "-", "#", "##", "//05/1987",
                                         "24/X/2013", "2015") | lych$identifiedby == 0] <- NA
lych <- lych %>% filter(!is.na(identifiedby))

#Correcting states' names
province <- lych$stateprovince
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
lych$stateprovince <- province
rm(province, lookup_states, get_states, i, j)
lych$stateprovince[lych$stateprovince == "?" | lych$stateprovince == "-"] <- NA
#plyr::count(lych$stateprovince) #Checking if everything has gone fine

#Filtering by states in which the campos rupestres occur 
#according to Silveira et al (2016) (19,624). NA's are not considered.
lych <- lych %>% filter(is.na(stateprovince) | stateprovince %in% c("Goias", 
                                                                  "Minas Gerais", 
                                                                  "Bahia", 
                                                                  "Pernambuco", 
                                                                  "Paraiba", 
                                                                  "Mato Grosso",
                                                                  "Distrito Federal"))

#Removing records without species lelych identification (19,014)
lych <- lych %>% filter(!is.na(species))
lych <- lych %>% filter(!species %in% c("sp.", "sp1"))
#plyr::count(lych$species)

#Removing records not based on preserved specimens and removing the 
#attribute 'basisofrecord' (18,983)
#plyr::count(lych$basisofrecord)
lych <- lych %>% filter(basisofrecord %in% c("Exsic", "PRESERVED_SPECIMEN",
                                           "s", "S", "PreservedSpecimen"))
lych <- lych %>% select(-basisofrecord)

#==============================================================================#

#=================================#
# CLEANING BY IDENTIFIER NAMES #
#=================================#

#Creating a vector to work with
identifiedby <- lych$identifiedby

#Generating a list of identifiers ordered by frequency of identifications (more identifications
#is supposed to be related with a better identification)
#Counting check
#names.count <- as.data.frame(plyr::count(identifiedby))
#names.count <- names.count[order(-names.count$freq), ]

#J. Semir
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("J.", "Semir"),
                              replace.by = "J. Semir",
                              not.replace = c("Seibert, J.B."),
                              replace = c("Semir, J; Loeuille, BFP", "Semir (UEC) 1992-11",
                                          "João, Semir", "Senir", "João Semir; b.1937; Semir",
                                          "Marques, D.; Semir, J.","Monge, M.; Semir, J.",
                                          "M.Monge & J.Semir"))

#B. Loueille
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("B.", "Loueille"),
                              replace.by = "B. Loueille",
                              not.replace = c("Lowille, B.I."),
                              replace = c("Benoit Loeuille","Luvielle","Loville",
                                          "Leville","Benoit, L","Benoit Loeuille (por foto)",
                                          "Benoit Loenille"))

#J. N. Nakajima
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("J. N.", "Nakajima"),
                              replace.by = "J. N. Nakajima",
                              not.replace = c(),
                              replace = c("J.Nakajima et R.Borges", "R. Borges & J. Nakajima",
                                          "Jimi Nakajima", "N, Nakajima", "Nakijima",
                                          "J. Nakajima et R.Borges", "J.N. Nakajima,R. Borges & M.M. Saavedra",
                                          "Jimi NakaJima", "R. Borges & J. Nakajima & M. Magenta"))

#C. M. Siniscalchi
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("C. M.", "Siniscalchi"),
                              replace.by = "C. M. Siniscalchi",
                              not.replace = c("C. M. Vieira"),
                              replace = c("fide Siniscalchi et al.", "G. M. Antar & C. M. Siniscalchi",
                                          "C. Siniscalchi (por foto)", "Carol Siniscalchi (por foto)"))

#D. J. N. Hindi
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("D. J. N.", "Hindi"),
                              replace.by = "D. J. N. Hindi",
                              not.replace = c(),
                              replace = c("Hind, DJN; Bautista, HP","D.J.N. Hind, September/",
                                          "D.J.N. Hind 1992", "D.J.N. Hind, July/","D.J.N. Hins, March/",
                                          "D.J.N. Hind, October/","N. Hind","D.J.N. Hind, March/",
                                          "D.J.N. Hind, September","D.J.N.Hind & H.P.Bautista",
                                          "Hind, D.J.N. & Bautista, H.P.","N. Hind & P. Bautista",
                                          "N. Hind (K) 2020-03-20","N. Hind (K) 2020-03-20"))

#N. F. F. MacLeish
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("N. F. F.", "MacLeish"),
                              replace.by = "N. F. F. MacLeish",
                              not.replace = c(),
                              replace = c("Nanda F. F. MacLeich","Nanda F. F. Macleish",
                                          "N.Mackeish","Nanda F. A Macleigh","N.F. MacLeish (GA)",
                                          "N.F. MacLeish (GA)","Nanda F. Fleming MacLeish; b.1953; MacLeish",
                                          "Nanda F.F. Macleish","Nanda F.F.Macleish","Nanda F, F, Macleish",
                                          "Nanda F.F.MacLiesh","Nanda F. F. MacLaish","Mackush, N.",
                                          "Nanda F. Fleming MacLeish"))

#H. Robinson
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("H.", "Robinson"),
                              replace.by = c("H. Robinson"),
                              not.replace = c("Robinson, BL","Robinson, E.A."),
                              replace = c("Harold Robinson em July.",
                                          "fide Robinson & Funk",
                                          "Harold Robinson em June."))

#N. Roque
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("N.", "Roque"),
                              replace.by = "N. Roque",
                              not.replace = c(),
                              replace = c("N.Roque & R.M.Liro","Roque, N; Oliveira, EC de"))

#G. M. Barroso
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("G. M.", "Barroso"),
                              replace.by = "G. M. Barroso",
                              not.replace = c("Barreto, M.","Barreto, M","G. M. Antar" ,
                                              "G. M. Barroso/ H. P. Bautista","Barros, R",
                                              "Barros, R."),
                              replace = c("Graziela Barroso","Graziela m. Barroso",
                                          "Graziela M. Barroso","G.M. Barroso & R. Marquete",
                                          "G.M.B","Graziela Maciel Barroso","Graziela M.Barroso"))

#J. B. Bringel
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("J. B.", "Bringel"),
                              replace.by = "J. B. Bringel",
                              not.replace = c(),
                              replace = c("João Bernardo Bringel Jr.","J.B.Bringel e D.A.Chaves",
                                          "Bringel; Nakajima, JA","Bringel Júnior, JBA; Chaves, DA"))

#R. A. X. Borges
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("R. A. X.", "Borges"),
                              replace.by = "R. A. X. Borges",
                              not.replace = c(),
                              replace = c("R.Borges & A.M.Teles","R. Borges & F. M. Ferreira"))

#M. Dematteis
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("M.", "Dematteis"),
                              replace.by = "M. Dematteis",
                              not.replace = c(),
                              replace = c("Almeida, G.; Dematteis, M."))

#D. A. Chaves
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.70, 
                              check.by = generate.names("D. A.", "Chaves"),
                              replace.by = "D. A. Chaves",
                              not.replace = c("Chaves, E.","Chaves, E","Chaves, ADM"),
                              replace = c())

#M. L. Guedes
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.70, 
                              check.by = generate.names("M. L.", "Guedes"),
                              replace.by = "M. L. Guedes",
                              not.replace = c("M. L. Gavilanes"),
                              replace = c("Guedes, ML; Passos Júnior, LA","Guedes, ML; Queiroz, EP",
                                          "Guedes, ML; Texeira, SR","Guedes, ML; Mattos Andrade, PE de",
                                          "Guedes, ML; França, TC","Guedes, ML; Prates, ARS",
                                          "Guedes, ML; Nascimento, AFS","Guedes, ML; Accioly, JAB",
                                          "Guedes, ML; Bautista, HP; Pinto, GCP","Guedes, ML; Noblick, LR; Faria, LSS; Paranhos, LH"))

#K. Calago
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.70, 
                              check.by = generate.names("K.", "Calago"),
                              replace.by = "K. Calago",
                              not.replace = c(),
                              replace = c())

#M. Monge
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M.", "Monge"),
                              replace.by = "M. Monge",
                              not.replace = c("M. M. Coelho"),
                              replace = c())
#H. F. Leitão-Filho
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("H. F.", "Leitão-Filho"),
                              replace.by = "H. F. Leitão-Filho",
                              not.replace = c(),
                              replace = c("H.F.L.Filho","Leitão Filho, HF; Semir, J",
                                          "H.F.L. FIlho","H.F. L. Filho","H.L.F. Filho",
                                          "Leitão f.; H.F.","H.F.L. Filho","H.F.L.. Filho"))

#G. Sancho
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("G.", "Sancho"),
                              replace.by = "G. Sancho",
                              not.replace = c(),
                              replace = c())

#H. A. Ogasawara
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("H. A.", "Ogasawara"),
                              replace.by = "H. A. Ogasawara",
                              not.replace = c(),
                              replace = c())

#R. S. Bianchini
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("R. S.", "Bianchini"),
                              replace.by = "R. S. Bianchini",
                              not.replace = c(),
                              replace = c("R. Simão-Bianchini","Simão Bianchini, R",
                                          "Simão Bianchini, R; Pastore, JA"))

#D. Marques
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("D.", "Marques"),
                              replace.by = "D. Marques",
                              not.replace = c("Marquete, R.","Marques, J.S."),
                              replace = c())

#J. E. Q. Faria
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("J. E. Q.", "Faria"),
                              replace.by = "J. E. Q. Faria",
                              not.replace = c(),
                              replace = c())

#R. C. Pereira
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("R. C.", "Pereira"),
                              replace.by = "R. C. Pereira",
                              not.replace = c("Pereira, BAS","E. Pereira",
                                              "Pereira, AMS","Pereira, LA",
                                              "Pereira, L.A.","L.C.Pereira",
                                              "G. C. Pereira Pinto"),
                              replace = c("R. Pereira & M.F. Cavalcanti","Rita Pereira; Cavalcanti, M.F." ,
                                          "\"\"Pereira, R.; Cavalcante, M.F.\"\"","Rita Pereira",
                                          "R. Pereira; M.F. Cavalcanti", "Rita Pereira, R"))

#R. Barros
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("R.", "Barros"),
                              replace.by = "R. Barros",
                              not.replace = c("Barbosa, E", "Barbosa, F.","Barreto, M",
                                              "Barreto, M."),
                              replace = c())

#J. B. Cândido
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("J. B.", "Cândido"),
                              replace.by = "J. B. Cândido",
                              not.replace = c("Cândido, E.S."),
                              replace = c())

#J. F. Pruski
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("J. F.", "Pruski"),
                              replace.by = "J. F. Pruski",
                              not.replace = c(),
                              replace = c("J. Pruski (MO)","John Pruski"))

#L. K. Kirkman
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("L. K.", "Kirkman"),
                              replace.by = "L. K. Kirkman",
                              not.replace = c(),
                              replace = c("L.Katherine Kirkman","L.. Katherine Kirkman"))

#R. L. Esteves
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("R. L.", "Esteves"),
                              replace.by = "R. L. Esteves",
                              not.replace = c(),
                              replace = c("Roberto L. Esteves"))

#V. C. Souza
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("V. C.", "Souza"),
                              replace.by = "V. C. Souza",
                              not.replace = c("Souza, F.S.","Souza, H.C. de","Souza, D.P.",
                                              "Souza-Silva, R.F.","Souza, H.F.","Souza, F.O. de",
                                              "Souza, FO"),
                              replace = c())

#B. M. T. Walter
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.65, 
                              check.by = generate.names("B. M. T.", "Walter"),
                              replace.by = "B. M. T. Walter",
                              not.replace = c(),
                              replace = c("Walter, BMT; Bringel Júnior, JBA","N. G. Ralha & B. M. T Walter"))

#E. K. O. Hattori
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("E. K. O.", "Hattori"),
                              replace.by = "E. K. O. Hattori",
                              not.replace = c(),
                              replace = c())

#C. E. B. Proença
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("C. E. B.", "Proença"),
                              replace.by = "C. E. B. Proença",
                              not.replace = c(),
                              replace = c())

#M. E. Mansanares
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M. E.", "Mansanares"),
                              replace.by = "M. E. Mansanares",
                              not.replace = c("Mansur, T." ),
                              replace = c())

#F. A. Santana
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("F. A.", "Santana"),
                              replace.by = generate.names("F. A.", "Santana")[1],
                              not.replace = c(),
                              replace = c())

#M. G. Staudt
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M. G.", "Staudt"),
                              replace.by = generate.names("M. G.", "Staudt")[1],
                              not.replace = c(),
                              replace = c())

#M. Barreto
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M.", "Barreto"),
                              replace.by = generate.names("M.", "Barreto")[1],
                              not.replace = c(),
                              replace = c())

#R. A. Pacheco
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("R. A.", "Pacheco"),
                              replace.by = generate.names("R. A.", "Pacheco")[1],
                              not.replace = c("Pacheco, LM"),
                              replace = c())

#J. G. Baker
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("J. G.", "Baker"),
                              replace.by = generate.names("J. G.", "Baker")[1],
                              not.replace = c("Baker, G.S."),
                              replace = c())

#L. Moura
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("L.", "Moura"),
                              replace.by = generate.names("L.", "Moura")[1],
                              not.replace = c(),
                              replace = c("Moura, L; Ogasawara, HA"))

#M. S. Oliveira
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M. S.", "Oliveira"),
                              replace.by = generate.names("M. S.", "Oliveira")[1],
                              not.replace = c("Oliveira, E.C.","Oliveira, EVS",
                                              "Oliveira, EC", "E. C. Oliveira",
                                              "Oliveira, CT","Oliveira, N.V.",
                                              "C. T. Oliveira","J. A. Oliveira",
                                              "Oliveira, RC","Oliveira, D.G.",
                                              "Oliveira, PA","Oliveira, EC de",
                                              "Oliveira, E.L.P.G. de","Oliveira, D." ,
                                              "Oliveira, CC","Oliveira, M.A."),
                              replace = c())

#A. M. Teles
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("A. M.", "Teles"),
                              replace.by = generate.names("A. M.", "Teles")[1],
                              not.replace = c(),
                              replace = c("Teles, AM; Nakajima, JN"))

#F. S. Freitas
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("F. S.", "Freitas"),
                              replace.by = generate.names("F. S.", "Freitas")[1],
                              not.replace = c("Freira, T.C."),
                              replace = c())

#C. Jeffrey
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("C.", "Jeffrey"),
                              replace.by = generate.names("C.", "Jeffrey")[1],
                              not.replace = c(),
                              replace = c("C.Jerrfry" ))

#A. Sciamarelli
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("A.", "Sciamarelli"),
                              replace.by = generate.names("A.", "Sciamarelli")[1],
                              not.replace = c(),
                              replace = c())

#H. P. Bautista
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("H. P.", "Bautista"),
                              replace.by = generate.names("H. P.", "Bautista")[1],
                              not.replace = c("Bautista, G."),
                              replace = c())

#J. Badini
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("J.", "Badini"),
                              replace.by = generate.names("J.", "Badini")[1],
                              not.replace = c(),
                              replace = c())

#R. C. Forzza
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("R. C.", "Forzza"),
                              replace.by = generate.names("R. C.", "Forzza")[1],
                              not.replace = c("R. C. Pereira"),
                              replace = c())

#M. F. Cavalcanti
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M. F.", "Cavalcanti"),
                              replace.by = generate.names("M. F.", "Cavalcanti")[1],
                              not.replace = c(),
                              replace = c())

#G. Hatschbach
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("G.", "Hatschbach"),
                              replace.by = generate.names("G.", "Hatschbach")[1],
                              not.replace = c(),
                              replace = c())

#S. Nunes
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("S.", "Nunes"),
                              replace.by = generate.names("S.", "Nunes")[1],
                              not.replace = c("Nunes, Y.R.F."),
                              replace = c())

#H. Schumacher
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("H.", "Schumacher"),
                              replace.by = generate.names("H.", "Schumacher")[1],
                              not.replace = c(),
                              replace = c())

#M. Alves
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("M.", "Alves"),
                              replace.by = generate.names("M.", "Alves")[1],
                              not.replace = c(),
                              replace = c())

#E. C. Oliveira
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("E. C.", "Oliveira"),
                              replace.by = generate.names("E. C.", "Oliveira")[1],
                              not.replace = c("Oliveira, CT","M. S. Oliveira",
                                              "Oliveira, PA","Oliveira, E.L.P.G. de",
                                              "Oliveira, EVS","Oliveira, N.V.",
                                              "J. A. Oliveira","Oliveira, D.","Oliveira, M.A.",
                                              "C. T. Oliveira","C.T. Oliveira","Oliveira, D.G.",
                                              "Oliveira, E.L.P.G. de no. 292","Oliveira, CC" ),
                              replace = c())

#J. A. Siqueira-Filho
identifiedby <- replace.names(x = identifiedby, top = 0.85, bottom = 0.7, 
                              check.by = generate.names("J. A.", "Siqueira-Filho"),
                              replace.by = generate.names("J. A.", "Siqueira-Filho")[1],
                              not.replace = c("Siqueira, L.C. de","J.C. Siqueira",
                                              "J. A. Oliveira"),
                              replace = c())

#names.count <- as.data.frame(plyr::count(identifiedby))
#names.count <- names.count[order(-names.count$freq), ]

#write.csv(names.count, file = "identificadores_lych.csv", row.names = FALSE)

#Replacing column
lych$identifiedby <- identifiedby

#Filtering (14,659)
specialists <- c("B. Loueille",
                 "J. Semir",
                 "D. J. N. Hindi",
                 "J. N. Nakajima",
                 "N. F. F. MacLeish",
                 "H. Robinson",
                 "G. M. Barroso",
                 "C. M. Siniscalchi",
                 "N. Roque",
                 "J. B. Bringel",
                 "H. F. Leitão-Filho",
                 "M. L. Guedes",
                 "D. A. Chaves",
                 "M. Dematteis",
                 "R. L. Esteves",
                 "R. C. Pereira",
                 "H. A. Ogasawara",
                 "J. F. Pruski",
                 "M. Monge",
                 "G. Sancho",
                 "J. B. Cândido",
                 "H. P. Bautista",
                 "A. M. Teles",
                 "D. Marques",
                 "J. E. Q. Faria",
                 "R. S. Bianchini",
                 "C. E. B. Proença",
                 "R. A. X. Borges",
                 "M. E. Mansanares",
                 "M. G. Staudt",
                 "E. K. O. Hattori",
                 "C. Jeffrey",
                 "R. A. Pacheco",
                 "L. Moura",
                 "F. S. Freitas",
                 "J. G. Baker",
                 "F. A. Santana",
                "M. Alves",
                 "E. C. Oliveira",
                 "R. C. Forzza",
                 "J. A. Siqueira-Filho")
lych <- lych %>% filter(identifiedby %in% specialists)

rm(identifiedby, specialists)

#====================================================================================================#

#===========================#
# CLEANING SCIENTIFIC NAMES #
#===========================#

#Generating a column for scientific names (without authors and including infraspecific epithet)
lych$gen_sp <- paste(lych$genus,
                    lych$species,
                    lych$subspecies,
                    sep = " ")

#taxa <- plyr::count(lych$gen_sp)
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
#write.csv(taxa_suggested, file = "lists/Lychnophorinae/taxa_suggested.csv", row.names = F)

#Loading *.csv after manual correction
taxa_corrected <- read.csv("lists/Lychnophorinae/taxa_corrected.csv", stringsAsFactors = F, 
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

#Removing invalid taxa (14,469)
lych$gen_sp <- gsub("NA", "", lych$gen_sp)
lych$gen_sp <- trimws(lych$gen_sp)
invalid_taxa <- taxa_gensp$replace[is.na(taxa_gensp$gen)]
lych <- lych[!lych$gen_sp %in% invalid_taxa, ]

#Correcting the dataset
taxa_gensp$replace <- taxa_corrected$taxa
taxa_gensp <- taxa_gensp[!is.na(taxa_gensp$gen), ]
for(i in 1:nrow(lych)){
  for(j in 1:nrow(taxa_gensp)){
    if(lych$gen_sp[i] == taxa_gensp$replace[j]){
      lych$gen_sp[i] <- paste(taxa_gensp$gen[j], taxa_gensp$sp[j], sep = " ")
    }
  }
}

#Checking if there is any last corrections to be made
lych$gen_sp[!lych$gen_sp %in% paste(taxa_gensp$gen, taxa_gensp$sp, sep = " ")]

rm(taxa_corrected, taxa_gensp, invalid_taxa, str)

#======================================================================================#

#================#
# SPLITTING DATA #
#================#

#Coords (6,667)
lych_coord <- lych %>% filter(!is.na(latitude) & !is.na(longitude))

#No coords (7,802)
lych_noCoord <- lych %>% filter(is.na(latitude) | is.na(longitude))

#======================================================================================#

#======================#
# CLEANING COORDINATES #
#======================#

#Cleaning (6,323)
lych_coordFlagged <- lych_coord %>% clean_coordinates(lon = "longitude",
                                                    lat = "latitude",
                                                    species = "gen_sp",
                                                    value = "flagged",
                                                    tests = c("equal", "gbif", 
                                                              "institutions", 
                                                              "outliers", "seas",
                                                              "zeros"))

invalid_coords <- lych_coord[lych_coordFlagged == FALSE, ]
lych_coordClean <- lych_coord[lych_coordFlagged  == TRUE, ]

#Binding invalid_coords to lych_noCoord (8,146)
lych_noCoord <- rbind(lych_noCoord, invalid_coords)

rm(invalid_coords, lych_coordFlagged, lych_coord)

#Writing *.csv 
#write.csv(lych_coordClean, file = "Datasets/lychloziaceae/lych_coordClean.csv", row.names = FALSE)
#======================================================================================#

#===============================#
# CLEANING MUNICIPALITIES NAMES #
#===============================#

#Dividing dataset in registers with and without values for stateprovince (1,130 and 7,016)
lych_noState <- lych_noCoord %>% filter(is.na(stateprovince))
lych_state <- lych_noCoord %>% filter(!is.na(stateprovince))

#NA values for municipalities (2,380)
lych_na <- lych_noCoord %>% filter(is.na(municipality_gbif) & is.na(county))

#Counting check
#names.count <- as.data.frame(plyr::count(lych_noState$county))
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
lych_noState$municipality_gbif_std <- correct.mun(lych_noState$municipality_gbif)
lych_state$municipality_gbif_std <- correct.mun(lych_state$municipality_gbif)
lych_noState$county_std <- correct.mun(lych_noState$county)
lych_state$county_std <- correct.mun(lych_state$county)

#Selecting mun_concla municipalities in which the campos rupestres occur
mun_concla_cr <- mun_concla[mun_concla$Nome_Município %in% list_mun_std, ]

#nrow(mun_concla_cr) is not equal to length(list_mun_std). Why?
#Municipalities from the São Paulo state
#list_mun_std[!list_mun_std %in% mun_concla_cr$Nome_Município]

#Manually checking and correcting by crossing the list assigned as 'check' with the mun_concla_cr dataset. 
#check <- data.frame(plyr::count(lych_noState$municipality_gbif_std)) #ok
lych_noState$municipality_gbif_std[lych_noState$municipality_gbif_std == 
                                    "Anaje"] <- "Anage"
lych_noState$municipality_gbif_std[lych_noState$municipality_gbif_std == 
                                    "Gouvea"] <- "Gouveia"
lych_noState$municipality_gbif_std[lych_noState$municipality_gbif_std == 
                                    "Sao Joao D'Alianca"] <- "Sao Joao D'alianca"
lych_noState$municipality_gbif_std[lych_noState$municipality_gbif_std 
                                  %in% c("Brasil", "Goias",
                                         "Mun?")] <- NA

#check <- data.frame(plyr::count(lych_noState$county_std)) #ok
lych_noState$county_std[lych_noState$county_std 
                       %in% c("Brasil", "5Km Sul da Cidade")] <- NA
lych_na <- rbind(lych_na, 
                lych_noState[which(is.na(lych_noState$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                    is.na(lych_noState$county_std)), 1:18]) 
lych_noState <- lych_noState[-which(is.na(lych_noState$municipality_gbif_std) & #removing registers with NA values for these two attributes from lych_noState
                                    is.na(lych_noState$county_std)), ]

#check <- data.frame(plyr::count(lych_state$municipality_gbif_std))
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std %in% 
                                   c("Alto Garca")] <- "Alto Garcas"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std %in% 
                                  c("Alto Paraiso","Alto Paraiso Goias")] <- "Alto Paraiso de Goias"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Abaira   Piata"] <- "Abaira"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Anaje"] <- "Anage"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Barra de Mendes"] <- "Barra do Mendes"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Betania/Floresta"] <- "Betania"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std %in% 
                                   c("Braslilia")] <- "Brasilia"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Catugi"] <- "Catuji"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Cajazeiras Sao Jose das Piranhas"] <- "Cajazeiras"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Cidade Ecletica"] <- "Santo Antonio do Descoberto"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Conceicao do Mato de Dentro"] <- "Conceicao do Mato Dentro"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Corumba  de Goias"] <- "Corumba de Goias"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std %in% 
                                   c("Delfinopolis (?)","Delfinopoilis")] <- "Delfinopolis"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Felixandia"] <- "Felixlandia"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Gouvea"] <- "Gouveia"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Olhos D#?#Agua"] <- "Olhos D'agua"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Planaltina de Goias"] <- "Planaltina"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Rio das Contas"] <- "Rio de Contas"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Santana do Pirapama"] <- "Santana de Pirapama"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Sao Tome das Letras"] <- "Sao Thome das Letras"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Sao Joao da Alianca"] <- "Sao Joao D'alianca"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std 
                                %in% c("Sao Joao Del'rei", "Sao Jose D'el Rei")] <- "Sao Joao Del Rei"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std == 
                                  "Terezina de Goias"] <- "Teresina de Goias"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std 
                                %in% c("Varzea de Palma")] <- "Varzea da Palma"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std 
                                %in% c("Vila Bela Santissima Trindade", "Vila Bela da Sma Trindade",
                                       "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
lych_state$municipality_gbif_std[lych_state$municipality_gbif_std 
                                %in% c("Chapada da Cotagem", "Chapada Diamantina",
                                       "Chapada do Araripe", "Chapada dos Guimaraes",
                                       "Chapada dos Veadeiros", "Paraiba", "Pernambuco",
                                       "Serra do Cipo", "Serra do Espinhaco", "Serra do Tombador",
                                       "Ufba Ondina", "Veadeiros", "Bahia", "Chapadao do Ceu", "Goias", 
                                       "Minas Gerais", "Serra do Cipo", "Serra do Espinhaco", 
                                       "Serra do Tombador", "Veadeiros", "Bahia", "Betim/Brumadinho",
                                       "Coletada A Margem da Estrada de Macambinha","Br 135 Km 404")] <- NA


#check <- data.frame(plyr::count(lych_state$county_std))
lych_state$county_std[lych_state$county_std %in% 
                       c("Alto Paraiso", "Alto Paraa-so de Goias")] <- "Alto Paraiso de Goias"
lych_state$county_std[lych_state$county_std == 
                       "Anaje"] <- "Anage"
lych_state$county_std[lych_state$county_std == 
                       "Barra de Mendes"] <- "Barra do Mendes"
lych_state$county_std[lych_state$county_std == 
                       "Br 020, 5Km Ne de Sobradinho"] <- "Sobradinho"
lych_state$county_std[lych_state$county_std == 
                       "Cristalina Mun"] <- "Cristalina"
lych_state$county_std[lych_state$county_std %in% 
                       c("Grao Mogol Mun")] <- "Grao Mogol"
lych_state$county_std[lych_state$county_std == 
                       "Joaquim Fela-cio"] <- "Joaquim Felicio"
lych_state$county_std[lych_state$county_std %in% 
                       c("Morro do Chapeu Mun")] <- "Morro do Chapeu"
lych_state$county_std[lych_state$county_std == 
                       "Rio das Contas"] <- "Rio de Contas"
lych_state$county_std[lych_state$county_std %in% 
                       c("Sao Joao da Alianca", 
                         "Sao Joao D##Alianca")] <- "Sao Joao D'alianca"
lych_state$county_std[lych_state$county_std %in% 
                       c("Sao Tome das Letras ")] <- "Sao Thome das Letras "
lych_state$county_std[lych_state$county_std %in% 
                       c("Vila Bela de Santa-ssima Trindade", 
                         "Vila Bela de Santissima Trindade",
                         "Vl Bela da Sma Trindade")] <- "Vila Bela da Santissima Trindade"
lych_state$county_std[lych_state$county_std == 
                       "Santana de Pirapama"] <- "Santana do Pirapama"
lych_state$county_std[lych_state$county_std == 
                       "Brasa-lia"] <- "Brasilia"
lych_state$county_std[lych_state$county_std 
                     %in% c("", "Chapada dos Guimaraes",
                            "Chapada dos Veadeiros", "Chapada Gaucha",
                            "Goias", "Minas Gerais", "No Disponible",
                            "Pe", "Ufba Ondina", "Veadeiros", 
                            "Minas Gerais", "Mato Grosso", "Serra de Ibitipoca",
                            "Serra do Cabral", "Serra do Cipo", "Serra do Espinhaco")] <- NA
lych_na <- rbind(lych_na, 
                lych_state[which(is.na(lych_state$municipality_gbif_std) & #binding registers with NA values for these two attributes
                                  is.na(lych_state$county_std)), 1:18]) 
lych_state <- lych_state[-which(is.na(lych_state$municipality_gbif_std) & #removing registers with NA values for these two attributes from lych_noState
                                is.na(lych_state$county_std)), ]

#Filtering lych_noState and lych_state with mun_concla_cr
lych_noState_filt <- lych_noState %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                              county_std %in% mun_concla_cr$Nome_Município)
lych_state_filt <- lych_state %>% filter (municipality_gbif_std %in% mun_concla_cr$Nome_Município |
                                          county_std %in% mun_concla_cr$Nome_Município)

#===================================================================================================#

#================================#
# INFERRING MUNICIPALITIES NAMES #
#================================#

#Standardizing 'locality'
lych_na$locality_std <- correct.mun(lych_na$locality)

#Vector with municipalities names from mun_concla to be used with grepl
grepl_munc_concla_cr <- c()
for(i in 1:nrow(mun_concla_cr)){
  grepl_munc_concla_cr[i] <- paste0("\\b", mun_concla_cr$Nome_Município[i],
                                    "\\b")
}

#Vector with correspondent municipalities names
vec <- 1:length(lych_na$locality_std)
for(i in 1:length(grepl_munc_concla_cr)){
  for(j in 1:length(lych_na$locality_std)){
    if(grepl(pattern = grepl_munc_concla_cr[i], x = lych_na$locality_std[j])){
      vec[j] <- mun_concla_cr$Nome_Município[i]
    }
  }
}
vec[vec %in% 1:length(lych_na$locality_std)] <- NA

#New column with inferred municipality names
lych_na$municipality <- vec

#Removing registers with NA for 'municipality'
lych_na <- lych_na %>% filter(!is.na(municipality))

#Concatenating information on municipality into a unique column 
lych_noState_filt$municipality <- NA
for(i in 1:nrow(lych_noState)){
 if(!is.na(lych_noState_filt$municipality_gbif_std[i])){
  lych_noState_filt$municipality[i] <- lych_noState_filt$municipality_gbif_std[i] 
} else if(!is.na(lych_noState_filt$county_std[i])){
  lych_noState_filt$municipality[i] <- lych_noState_filt$county_std[i]
}
}

lych_state_filt$municipality <- NA
for(i in 1:nrow(lych_state)){
  if(!is.na(lych_state_filt$municipality_gbif_std[i])){
    lych_state_filt$municipality[i] <- lych_state_filt$municipality_gbif_std[i] 
  } else if(!is.na(lych_state_filt$county_std[i])){
    lych_state_filt$municipality[i] <- lych_state_filt$county_std[i]
  }
}

#Concatenating data sets
lych_noState_filt <- lych_noState_filt[ , -which(names(lych_noState_filt) 
%in% c("county","county_std",
      "municipality_gbif",
     "municipality_gbif_std"))]

lych_state_filt <- lych_state_filt[ , -which(names(lych_state_filt) 
                                           %in% c("county","county_std",
                                                  "municipality_gbif",
                                                  "municipality_gbif_std"))]
lych_na <- lych_na[ , -which(names(lych_na) 
                           %in% c("county", "municipality_gbif", "locality_std"))]

lych_noCoord_inf<- rbind(lych_noState_filt, lych_state_filt, lych_na)

#Registers occurring in homonyms municipalities
reg_hom <- lych_noCoord_inf[lych_noCoord_inf$municipality %in% homonyms, ]
reg_hom <- reg_hom %>% filter(!is.na(stateprovince))
lych_noCoord_inf <- lych_noCoord_inf %>% filter(!municipality %in% homonyms)
lych_noCoord_inf <- rbind(lych_noCoord_inf, reg_hom) #binding them after removing those without state/province information

rm(lych_na, lych_noState, lych_noState_filt, lych_state, lych_state_filt,
   mun_concla, grepl_munc_concla_cr, list_mun, list_mun_std, vec, cr,
   mun, crswgs84, lych_noCoord, reg_hom)

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
for(i in 1:length(lych_noCoord_inf$municipality)){
  for(j in 1:length(mun_centr_cr$`minor area`)){
    if(lych_noCoord_inf$municipality[i] == mun_centr_cr$`minor area`[j]){
      lych_noCoord_inf$latitude[i] <- mun_centr_cr$`centroid (lat)`[j]
      lych_noCoord_inf$longitude[i] <- mun_centr_cr$`centroid (long)`[j] 
    } 
  }
}

#Homonyms
for(i in 1:length(lych_noCoord_inf$municipality)){
  for(j in 1:length(centr_hom$`minor area`)){
    if(lych_noCoord_inf$municipality[i] == centr_hom$`minor area`[j] &
       lych_noCoord_inf$stateprovince[i] == centr_hom$`major area`[j] &
       !is.na(lych_noCoord_inf$stateprovince[i])){
      lych_noCoord_inf$latitude[i] <- centr_hom$`centroid (lat)`[j]
      lych_noCoord_inf$longitude[i] <- centr_hom$`centroid (long)`[j] 
    } 
  }
}

#Removing NA's for coordinates
lych_noCoord_inf <- lych_noCoord_inf %>% filter(!is.na(latitude) & !is.na(longitude))

rm(mun_centr, mun_centr_br, mun_centr_cr, mun_concla_cr, centr_hom, homonyms)

#=================================================================================================#

#================#
# CLEAN DATA SET #
#================#

#Concatenating data sets and removing unnimportant columns
lych_noCoord_inf <- lych_noCoord_inf %>% select(id,
                                              institutioncode,
                                              collectioncode,
                                              catalognumber,
                                              gen_sp,
                                              subspecies,
                                              identifiedby,
                                              latitude,
                                              longitude)

lych_coordClean <- lych_coordClean %>% select(id,
                                            institutioncode,
                                            collectioncode,
                                            catalognumber,
                                            gen_sp,
                                            subspecies,
                                            identifiedby,
                                            latitude,
                                            longitude)

lych_clean <- rbind(lych_coordClean, lych_noCoord_inf)

#Is there any duplicated registers? 
#lych_clean[duplicated(lych_clean$id) == TRUE, ] #no

rm(lych, lych_coordClean, lych_noCoord_inf, i, j, correct.mun, generate.names,
   replace.names, titling)

#Writing *.csv 
#write.csv(lych_clean, file = "Datasets/Lychnophorinae/lych_clean.csv", row.names = FALSE)

#====================================================================================================#

#======================#
# DATASET - CR RECORDS #
#======================#

#Reading clean dataset
lych_clean <- read.csv(file = "datasets/Lychnophorinae/lych_clean.csv", na.strings = c("", NA))

#Standardizing gen_sp column
lych_clean$gen_sp <- gsub(" ", "_", lych_clean$gen_sp)

#Reading shapefiles: campos rupestres and spatial grids
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
grids_cr <- readOGR("shapefiles/grids_cr/grids_cr.shp") 

#Projecting
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(cr) <- crswgs84

#Intersecting coordinates with the cr shapefile and, then, with the grids 
coords <- lych_clean 
coordinates(coords) <- ~ longitude + latitude
proj4string(coords) <- crswgs84
coords_2 <- over(coords, cr)
coords_2$id_2 <- 1:nrow(coords_2)
coords <- lych_clean
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
write.csv(coords, "datasets/Lychnophorinae/lych_cr.csv", row.names = FALSE)
