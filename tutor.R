install.packages("raster")
install.packages("rasterVis")
install.packages("maptools")
install.packages("rworldmap")
install.packages("googleMaps") ##no funciona para versión 3
install.packages("dismo")

library(sp)
library(raster)
library(rasterVis)
library(maptools)
library(rworldmap)
library(googleVis)
library(googleMaps)
library(dismo)

newmap <- getMap(resolution="coarse")
plot(newmap)
mapCountryData()
mapCountryData(mapRegion="europe")
mapGriddedData()
mapGriddedData(mapRegion="europe")

data(Exports)
View(Exports)
Geo <- gvisGeoMap(Exports,locationvar="Country",
                  numvar="Profit",
                  options=list(height=400,dataMode='regions'))

Geo <- gvisGeoMap(Exports, locationvar="Country", numvar="Profit", 
                  options=list(height=400, dataMode='regions'))

plot(Geo)
data(Andrew)
?Andrew

M1 <- gvisMap(Andrew,"LatLong","Tip",
              options=list(showTip=TRUE,showLine=FALSE,
                           enableScrollWheel=TRUE,
                           mapType='satellite',
                           useMapTypeControl=TRUE,
                           width=800,
                           height=400))

M1 <- gvisMap(Andrew, "LatLong" , "Tip",
              options=list(showTip=TRUE, showLine=F,
                           enableScrollWheel=TRUE,
                           mapType='satellite', 
                           useMapTypeControl=TRUE, 
                           width=800,height=400))
plot(M1)

Andrew

mymap <- gmap("México")
plot(mymap)
?gmap

mymap <- gmap("México",type="satellite")
plot(mymap)

mymap <- gmap("México",type="satellite",exp=3)
plot(mymap)

mymap <- gmap("Europe",exp=1)
plot(mymap)