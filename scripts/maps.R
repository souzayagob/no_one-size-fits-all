#Function to get equivalent color with different transparency. From: https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}

#Show colors
library(scales)
show_col(c("darkgreen","darkblue","darkred","darkorange"))
show_col(c("#231151FF",
           "#B63679FF",
           "#4AC16DFF",
           "#CFE11CFF"))

#Colours (NE, SE, BCP, SC)
show_col(c("#008805",
           "#0000A3",
           "#FF0005",
           "#FF7F00"))
#======================================================================================================#

library(rgdal)
library(raster)
library(viridis)
library(rgeos)
library(latticeExtra)

#Defining projection
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

#Reading shapefiles
cr <- readOGR("shapefiles/shapefilecamporupestre/cr.shp") #Silveira (2016)
br <- readOGR("shapefiles/br_unidades_da_federacao/BRUFE250GC_SIR.shp", encoding = "UTF-8") #IBGE: https://downloads.ibge.gov.br/downloads_geociencias.htm

#Projecting
proj4string(cr) <- crswgs84
br <- spTransform(br, crswgs84)

#br_campos_rup
br_cr <- raster(extent(-50.5, -38.5, -23.25, -8.5),
                       crs = crswgs84, nrows = 1, ncols = 1)
br_cr <- rasterToPolygons(br_cr)
plot(br)
plot(br_cr, col = "#BEBEBE7F", add = TRUE)
plot(cr, col = "black", add = T)

#Espinhaco
campos_rup <- crop(cr, extent(c(-50.5, -38.5, -23.25, -8.5)))
br_mg <- br[br$NM_ESTADO %in% c("MINAS GERAIS", "SÃO PAULO"), ]
br_ba <- br[br$NM_ESTADO %in% c("BAHIA"), ]
br_go <- br[br$NM_ESTADO %in% c("DISTRITO FEDERAL", "GOIÁS"), ]
campos_rup@data$sector <- c() 
for(i in 1:nrow(campos_rup)){
  if(gIntersects(campos_rup[i, ], br_mg)){
    campos_rup@data$sector[i] <- "Southern Espinhaço"
  } else if(gIntersects(campos_rup[i, ], br_ba)){
    campos_rup@data$sector[i] <- "Northern Espinhaço"
  } else if(gIntersects(campos_rup[i, ], br_go)){
    campos_rup@data$sector[i] <- "Brazilian Central Plateau"
  }
}

campos_rup$sector <- factor(campos_rup$sector, 
                           c("Southern Espinhaço", "Northern Espinhaço", 
                             "Brazilian Central Plateau"),
                           ordered = is.ordered(campos_rup$sector))

#elevation raster
alt <- raster('HYP_HR_SR/HYP_HR_SR.tif')
e <- as(extent(-50.5, -38.5, -23.25, -8.5), 'SpatialPolygons')
crs(e) <- crswgs84
r <- crop(alt, e)
rr <- mask(r, br)

campos_rup$color_in <- c()
for(i in 1:nrow(campos_rup@data)){
  if(campos_rup$sector[i] == "Southern Espinhaço"){
    campos_rup$color_in[i] <- makeTransparent("#FF0005", alpha = 0.65)
  } else if(campos_rup$sector[i] == "Northern Espinhaço"){
    campos_rup$color_in[i] <- makeTransparent("#008805", alpha = 0.65)
  } else if(campos_rup$sector[i] == "Brazilian Central Plateau"){
    campos_rup$color_in[i] <- makeTransparent("#FF0005", alpha = 0.65)
  }
}

campos_rup$sector <- factor(campos_rup$sector, levels = c("Southern Espinhaço",
                                                        "Northern Espinhaço",
                                                        "Brazilian Central Plateau"))

cr_plot  <- spplot(campos_rup, zcol = "sector",
       xlim = c(-50.5, -38.5), ylim = c(-23.25, -8.5), 
       sp.layout = list(list(br, "black", 
                             fill = "transparent", alpha = 0.7)),
       colorkey = FALSE,
       #colorkey = list(tck = 0),
       scales = list(draw = FALSE),
       #scales = list(draw = TRUE, x = list(at = c(-46, -42, -38)),
       #              y = list(at =c(-19, -15, -11))),
       col = campos_rup$color_in,
       col.regions = makeTransparent(c("#FF0005", "#008805", "#FF0005"), alpha = 0.5)) + as.layer(spplot(rr, col.regions = grey(1:100/100, alpha = 0.7), 
                                                                maxpixels = 2e10), 
                                                         under = TRUE)

#png("plots/cr_plot3.png",
    #height = 4, width = 4, units = 'in', res=300); cr_plot; dev.off()

#png("plots/br_plot.png",
    #height = 4, width = 4, units = 'in', res=300)
#plot(br)
#plot(br_cr, col = "#BEBEBE7F", add = TRUE)
#plot(cr, col = "black", add = TRUE)
#dev.off()
