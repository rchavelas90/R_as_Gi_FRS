##http://plotkml.r-forge.r-project.org/gpxbtour.html

library(XML)
readGPX <- function(
 gpx.file,
 metadata = TRUE,
 bounds = TRUE,
 waypoints = TRUE, 
 tracks = TRUE,  
 routes = TRUE   
)
 
{    
 options(warn = -1)    
 
 if(metadata==TRUE) { metadata <- .readGPX.element(gpx.file, "name") }    
 if(bounds==TRUE) { bounds <- .readGPX.element(gpx.file, "bounds") }    
 if(waypoints==TRUE) { waypoints <- .readGPX.element(gpx.file, "wpt") }
 if(tracks==TRUE) { tracks <- .readGPX.element(gpx.file, "trk") }
 if(routes==TRUE) { routes <- .readGPX.element(gpx.file, "rte") }
 
 gpx <- list(metadata=metadata, bounds=bounds, waypoints=waypoints, tracks=tracks, routes=routes)
 return(gpx)
}

## Read various elements from a *.gpx file:

.readGPX.element <- function(gpx.file, element) {
 # element = "metadata", "wpt", "rte", "trk"
 
 ret <- xmlTreeParse(gpx.file, useInternalNodes = TRUE)
 # top structure: 
 top <- xmlRoot(ret)
 
 # check if there is any content:
 if(any(grep(element, names(top)))) {
  
  # tracks:
  if(element=="trk"){   
   ret <- NULL
   nu <- which(names(top) %in% element)
   for(c in seq_along(nu)){
    lst <- which(names(top[[nu[c]]]) %in% "trkseg")
    nm <- names(top[[nu[c]]][[lst[1]]][[1]])
    ret[[c]] <- list(NULL)
    for(i in seq_along(lst)) {
     trkpt <- top[[nu[c]]][[lst[i]]]
     ret[[c]][[i]] <- data.frame(NULL)
     ## get columns (http://www.topografix.com/GPX/1/1/#type_wptType)
     lon <- as.numeric(xmlSApply(trkpt, xmlGetAttr, "lon"))
     lat <- as.numeric(xmlSApply(trkpt, xmlGetAttr, "lat"))
     ret[[c]][[i]][1:length(lon),"lon"] <- lon
     ret[[c]][[i]][1:length(lat),"lat"] <- lat
     if(!nm[[1]]=="NULL"){
      for(j in 1:length(nm)){
       xm <- as.character(sapply(sapply(xmlChildren(trkpt), function(x) x[[nm[[j]]]]), xmlValue))
       ret[[c]][[i]][1:length(xm), nm[[j]]] <- xm 
      }
     } 
    }
    names(ret[[c]]) <- xmlValue(top[[nu[c]]][["name"]])
   }   
  }
  
  if(element=="wpt"){
   ret <- data.frame(NULL)
   nu <- which(names(top) %in% element)
   nm <- names(top[[nu[1]]])
   for(i in seq_along(nu)) {
    # coordinates:
    ret[i, "lon"] <- as.numeric(xmlGetAttr(top[[nu[i]]], "lon"))
    ret[i, "lat"] <- as.numeric(xmlGetAttr(top[[nu[i]]], "lat"))
    if(!nm[[1]]=="NULL"){
     for(j in 1:length(nm)){
      ret[i, nm[[j]]] <- xmlValue(xmlChildren(top[[nu[i]]])[[nm[[j]]]])
     }  
    }
   }
  }
  
  if(element=="rte"){
   ret <- NULL
   nu <- which(names(top) %in% element)
   for(c in seq_along(nu)){
    ret[[c]] <- data.frame(NULL)
    lst <- which(names(top[[nu[c]]]) %in% "rtept")
    nm <- names(top[[nu[c]]][[lst[1]]])
    for(i in seq_along(lst)) {
     rtept <- top[[nu[c]]][[lst[i]]]
     ret[[c]][i, "lon"] <- as.numeric(xmlGetAttr(rtept, "lon"))
     ret[[c]][i, "lat"] <- as.numeric(xmlGetAttr(rtept, "lat"))
     if(!nm[[1]]=="NULL"){
      for(j in c("name","cmt","desc","sym","type")){
       try(ret[[c]][i, j] <- xmlValue(rtept[[j]]), silent = TRUE)
      }
     } 
    }
    names(ret)[c] <- xmlValue(top[[nu[c]]][["name"]])
   }
  }
  
  # bounds
  if(element=="bounds"){
   nu <- which(names(top) %in% element)
   ret <- matrix(rep(NA, 4), nrow=2, dimnames = list(c("lat", "lon"), c("min", "max")))
   # coordinates:
   ret[1,1] <- as.numeric(xmlGetAttr(top[[nu[1]]], "minlon"))
   ret[1,2] <- as.numeric(xmlGetAttr(top[[nu[1]]], "maxlon"))
   ret[2,1] <- as.numeric(xmlGetAttr(top[[nu[1]]], "minlat"))
   ret[2,2] <- as.numeric(xmlGetAttr(top[[nu[1]]], "maxlat"))
  }
  
  # metadata
  if(element=="name"){
   lst <- c("name","desc","author","email","url","urlname","time")
   nu <- which(names(top) %in% lst)
   if(!nu[[1]]=="NULL"){      
    ret <- data.frame(NULL)
    for(i in seq_along(lst)) {
     try(ret[1,lst[i]] <- xmlValue(top[[nu[[i]]]]), silent = TRUE)
    }
   }
  }
  
 }
 else { ret <- NULL }
 
 return(ret)
}

# enf of script;
getwd()
data <- readGPX("2013-09-17-19-53.gpx",metadata=F,bounds=F,waypoints=F,route=F)
setwd("GPX")

str(data)
str(data$tracks)
data2 <- data$tracks[[1]][[1]]
head(data2)
head(data2$time)
data2$delta <- c()
data2$time <- as.POSIXct(data2$time, format="%Y-%m-%dT%H:%M:%SZ")
plot(data2$lat,data2$lon)
data2$loc=paste(data2$lat,data2$lon,sep=":")

library(RgoogleMaps)
MyMap <- MapBackground(lat=data2$lat,lon=data2$lon)
PlotOnStaticMap(MyMap=MyMap,lat=data2$lat,lon=data2$lon,FUN=lines,destfile="ejm.png")


library(googleVis)
g <- gvisGeoChart(data2,"loc",options=list(region="MX"))
plot(g)
gvisGe

head(data2)
g2 <- gvisMap(data2,"loc",options=list(showLine=TRUE))
plot(g2)