##http://cnr.lwlss.net/GarminR/
# Read in xml files from Garmin ForeRunner
library(XML)
setwd("GPX")
filename="activity_61755378"


doc=xmlParse(paste(filename,".gpx",sep=""),useInternalNodes=TRUE)

# We need to extract the simple tree structure containing 
# the data which we can pass to xmlToDataFrame
top=xmlRoot(doc)
title=toString.XMLNode(top[[2]][[1]][[1]])

description=toString.XMLNode(top[[2]][[2]][[1]])
data=toString.XMLNode(top[[2]][[3]])
# All these tests in case description node is empty!
# When it is, data node is moved up one
# Really need to access nodes by name, but can't see how...
if(data=="NULL") data=toString.XMLNode(top[[2]][[2]])
if(data=="NULL") data=toString.XMLNode(top[[2]][[1]])

# Fill a data frame with interesting data
df=as.data.frame(xmlToDataFrame(data,c("numeric", "character", "integer")))
#Get the attributes of the XMLNode
if(toString.XMLNode(top[[2]][[3]])=="NULL") {
 if(toString.XMLNode(top[[2]][[2]])=="NULL") {
  attribs=xmlSApply(top[[2]][[1]],xmlAttrs)
 }else{
  attribs=xmlSApply(top[[2]][[2]],xmlAttrs)
 }
}else{
 attribs=xmlSApply(top[[2]][[3]],xmlAttrs)
}

df$lon=as.numeric(attribs[1,])
df$lat=as.numeric(attribs[2,])
colnames(df)=c("Elevation","DateTime","HeartRate","Longitude","Latitude")
str(df)
df$Elevation=as.numeric(df$Elevation)
df$HeartRate=as.integer(df$HeartRate)
df$DateTime=as.character(df$DateTime)
# Convert timestamp to number of seconds since start of run

date=substr(df$DateTime[1],1,10)
Time=substr(df$DateTime,12,19)
T0=strptime(Time[1],"%H:%M:%S")
Time=as.numeric(strptime(Time,"%H:%M:%S")-T0)
df$Seconds=Time

##R#Look for a pattern in the time of the records
rm(dfRic)
dfRic <- data.frame(Time=Time,Delta=0)
head(dfRic,10)
i <- 1
for(i in 1:(length(dfRic[[1]])-1)){
dfRic[i+1,2] <- dfRic[i+1,1]-dfRic[i,1]
head(dfRic,10)
}
head(dfRic,12)
hist(dfRic$Delta)
hist(dfRic$Delta)
plot(dfRic$Delta,type="l")


# Initialise columns
df$dNorth=0; df$dEast=0; df$dUp=0;
df$North=0; df$East=0; df$dDist=0; 
df$dDist2D=0; df$Dist2D=0

# Haversine formula is appropriate for calculating distances from lat/long
EarthRad=6371000
haverDist<-function(aLong,aLat,bLong,bLat){
 dLat=2*pi*(bLat-aLat)/360.0; dLon=2*pi*(bLong-aLong)/360.0
 a=(sin(dLat/2))^2+cos(2*pi*aLat/360)*cos(2*pi*bLat/360)*(sin(dLon/2)^2)
 return(EarthRad*2*atan2(sqrt(a),sqrt(1-a)))
}

# Calculate northings and eastings
df$East=haverDist(df[1,"Longitude"],df[1,"Latitude"],df$Longitude,df[1,"Latitude"])*sign(df$Longitude-df[1,"Longitude"])
df$North=haverDist(df[1,"Longitude"],df[1,"Latitude"],df[1,"Longitude"],df$Latitude)*sign(df$Latitude-df[1,"Latitude"])

# Calculate changes in position for each dt
for (x in 2:(length(df$DateTime)-1)) {
 sEast=sign(df[x,"Longitude"]-df[1,"Longitude"])
 sNorth=sign(df[x,"Latitude"]-df[1,"Latitude"])
 df$dEast[x]=sEast*haverDist(df[x-1,"Longitude"],df[1,"Latitude"],df[x,"Longitude"],df[1,"Latitude"])
 df$dNorth[x]=sNorth*haverDist(df[1,"Longitude"],df[x-1,"Latitude"],df[1,"Longitude"],df[x,"Latitude"])
 df$dUp[x]=df$Elevation[x]-df$Elevation[x-1]
 # 2D distance (ignoring hills)
 df$dDist2D[x]=haverDist(df[x-1,"Longitude"],df[x-1,"Latitude"],df[x,"Longitude"],df[x,"Latitude"])
}

df$dDist=sqrt(df$dNorth^2+df$dEast^2+df$dUp^2)
df$Dist=cumsum(df$dDist)
df$Dist2D=cumsum(df$dDist2D)

# Fit a spline function to the GPS coordinates & elevation
east=splinefun(df$Seconds,df$East)
north=splinefun(df$Seconds,df$North)
up=splinefun(df$Seconds,df$Elevation)
dist=splinefun(df$Seconds,df$Dist)
hr=approxfun(df$Seconds,df$HeartRate) # Some gaps in heart rate record, linear interpolation more robust

# Do finite centred differencing to give smoothest rate/gradient estimates
df$Speed=rep(0,length(df$Seconds))
df$Gradient=rep(0,length(df$Seconds))
for(x in 2:(length(df$Seconds)-1)){
 Dt=df[x+1,"Seconds"]-df[x-1,"Seconds"]
 Dd=df[x+1,"Dist"]-df[x-1,"Dist"]
 df[x,"Speed"]=Dd/Dt # m/s
 df[x,"Gradient"]=(df[x+1,"Elevation"]-df[x-1,"Elevation"])/Dd # m/m
}
df[1,"Speed"]=df[2,"Speed"]
df[length(df$Seconds),"Speed"]=df[length(df$Seconds)-1,"Speed"]
df[1,"Gradient"]=df[2,"Gradient"]
df[length(df$Seconds),"Gradient"]=df[length(df$Seconds)-1,"Gradient"]

# Smooth speed as it is unrealistically noisy
df$Speed=smooth(df$Speed)

# Fit a spline function to rate
speed=splinefun(df$Seconds,df$Speed)
pace<-function(t) sapply(1/speed(t),max,0)
ppace<-function(t) 1000*pace(t)/60

# Update dataframe with speed and pace
df$Speed=speed(df$Seconds)
df$Pace=pace(df$Seconds)

# Generate some plots
reportfile=paste(title,filename,".pdf",sep="")
print(paste("Building",reportfile))
pdf(reportfile)

# Generate time interpolation points
Num=2001
minT=0; maxT=max(df$Seconds)
interT=minT+(maxT-minT)*(0:Num)/Num

# Create a colour function for plots
colfunc=colorRampPalette(c("navy","white", "red3"),space="Lab")
cp=colfunc(500)
getCol<-function(colFrac) cp[1+round(499*colFrac)]

# Generate fractional variables for colouring plots
hrFrac=(hr(interT)-min(hr(interT)))/(max(hr(interT))-min(hr(interT)))
upFrac=(up(interT)-min(up(interT)))/(max(up(interT))-min(up(interT)))
pmax=min(c(60*7/1000,max(pace(interT)))) # Else scales ruined by stopping and walking
pFrac=(pace(interT)-min(pace(interT)))/(pmax-min(pace(interT)))

# Calculate Color Scales
hrLevels=min(hr(interT))+(1:length(cp))*(max(hr(interT))-min(hr(interT)))/length(cp)
upLevels=min(up(interT))+(1:length(cp))*(max(up(interT))-min(up(interT)))/length(cp)
pmax=min(c(7,max(ppace(interT)))) # Else scales ruined by stopping and walking
pLevels=min(ppace(interT))+(1:length(cp))*(pmax-min(ppace(interT)))/length(cp)

# Make a plotting dataframe, calculate displacement during each timestep
plt=data.frame(time=interT,east=east(interT),north=north(interT),up=up(interT),hr=hr(interT),distance=sapply(interT,dist),speed=speed(interT),pace=pace(interT))

# Elevation trace
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
plot(NULL,xlab="East (m)",ylab="North (m)",xlim=c(min(df$East),max(df$East)),ylim=c(min(df$North),max(df$North)),main=paste(title,date))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
points(east(interT),north(interT),pch=16,cex=0.6,col=getCol(upFrac))
# Draw legend
image(1, upLevels,matrix(data=upLevels, ncol=length(upLevels),nrow=1),col=cp,xlab="",ylab="Elevation (m)",xaxt="n")
layout(1)

# Heart rate trace
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
plot(NULL,xlab="East (m)",ylab="North (m)",xlim=c(min(df$East),max(df$East)),ylim=c(min(df$North),max(df$North)),main=paste(title,date))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
points(east(interT),north(interT),pch=16,cex=0.6,col=getCol(hrFrac))
# Draw legend
image(1, hrLevels,matrix(data=hrLevels, ncol=length(hrLevels),nrow=1),col=cp,xlab="",ylab="Heart Rate (bpm)",xaxt="n")
layout(1)

# Pace trace
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
plot(NULL,xlab="East (m)",ylab="North (m)",xlim=c(min(df$East),max(df$East)),ylim=c(min(df$North),max(df$North)),main=paste(title,date))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
points(east(interT),north(interT),pch=16,cex=0.6,col=getCol(pFrac))
# Draw legend
image(1, pLevels,matrix(data=pLevels, ncol=length(pLevels),nrow=1),col=cp,xlab="",ylab="Pace (min/km)",xaxt="n")
layout(1)

#####
op<-par(mfrow=c(2,2))

# Elevation timecourse
plot(df$Seconds/60,df$Elevation,xlab="Time (min)",ylab="Elevation (m)",type="l",lwd=2,col="red")

# Heart rate timecourse
plot(df$Seconds/60,df$HeartRate,xlab="Time (min)",ylab="Heart Rate (bpm)",type="l",lwd=2,col="red")

# Distance timecourse
#plot(df$Seconds/60,df$Dist/1000,xlab="Time (min)",ylab="Distance (km)",type="l",lwd=2,col="red")

# Speed timecourse
plot(df$Seconds/60,60*df$Speed/1000,xlab="Time (min)",ylab="Speed (km/min)",type="l",lwd=2,col="red")

# Pace timecourse
pmin=max(0,1000*min(df$Pace)/60)
pmax=min(7,1000*max(df$Pace)/60)
plot(df$Seconds/60,1000*df$Pace/60,xlab="Time (min)",ylab="Pace (min/km)",type="l",lwd=2,col="red",ylim=c(pmin,pmax))
title("Performance statistics with time (min)",line=-2,outer=TRUE)
par(op)

#####
op<-par(mfrow=c(2,2))

# Elevation timecourse
plot(df$Dist/1000,df$Elevation,xlab="Distance (km)",ylab="Elevation (m)",type="l",lwd=3,col="blue")

# Heart rate timecourse
plot(df$Dist/1000,df$HeartRate,xlab="Distance (km)",ylab="Heart Rate (bpm)",type="l",lwd=2,col="blue")

# Distance timecourse
#plot(df$Dist/1000,df$Dist/1000,xlab="Distance (km)",ylab="Distance (km)",type="l",lwd=2,col="blue")

# Speed timecourse
plot(df$Dist/1000,60*df$Speed/1000,xlab="Distance (km)",ylab="Speed (km/min)",type="l",lwd=2,col="blue")

# Pace timecourse
pmin=max(0,1000*min(df$Pace)/60)
pmax=min(7,1000*max(df$Pace)/60)
plot(df$Dist/1000,1000*df$Pace/60,xlab="Distance (km)",ylab="Pace (min/km)",type="l",lwd=2,col="blue",ylim=c(pmin,pmax))
title("Performance statistics with distance (km)",line=-2,outer=TRUE)
par(op)

#####
op<-par(mfrow=c(1,2))

hist(1000*plt$pace/60,breaks=21,xlab="Pace (min/km)",ylab="Frequency",main="")
hist(plt$hr,breaks=61,xlab="Heart Rate (bpm)",ylab="Frequency",main="")

title("Frequency histograms",line=-2,outer=TRUE)
par(op)

#####
op<-par(mfrow=c(2,2))

# Pace gradient correlation
gpCor=formatC(cor(df$Gradient,1000*df$Pace/60), digits=4)
plot(df$Gradient,1000*df$Pace/60,col="red",pch=16,xlab="Gradient",ylab="Pace (min/km)",main=paste("Correlation:",gpCor))

# Heart-rate gradient correlation
# need to strip out warming up period
minHR=mean(df$HeartRate)-1.98*sd(df$HeartRate)
times=df$Seconds[df$HeartRate>=minHR]
mint=min(times)
dfHR=df[df$Seconds>mint,]
gpCor=formatC(cor(dfHR$Gradient,dfHR$HeartRate), digits=4)
plot(dfHR$Gradient,dfHR$HeartRate,col="red",pch=16,xlab="Gradient",ylab="Heart Rate (bpm)",main=paste("Correlation:",gpCor))

# Pace time correlation
gpCor=formatC(cor(df$Seconds,1000*df$Pace/60), digits=4)
plot(df$Seconds,1000*df$Pace/60,col="red",pch=16,xlab="Time (s)",ylab="Pace (min/km)",main=paste("Correlation:",gpCor),
     ylim=c(1000*min(df$Pace)/60,min(c(7,1000*max(df$Pace)/60))))

# Heart-rate time correlation
gpCor=formatC(cor(dfHR$Seconds,dfHR$HeartRate), digits=4)
plot(dfHR$Seconds,dfHR$HeartRate,col="red",pch=16,xlab="Time (s)",ylab="Heart Rate (bpm)",main=paste("Correlation:",gpCor))

par(op)

dev.off()