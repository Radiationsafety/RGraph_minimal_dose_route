# This code is used for interpolation of measerements and for creation of a dose rate grid
library(automap); library(geosphere)  
library(raster); library(rasterVis);library(rgdal);library(geoR); library(gstat)

# Data file with measerements of dose date
inputdatafilename<-"data_RMO_OTZPVH_test.dat"
dataRMO<-read.table(inputdatafilename)  
colnames(dataRMO)=c("x","y","z") 

basemap<-raster("industrial site.gif") 
wsp<-extent(basemap)
# step of grid
xdist=distm (c(ymin(wsp), xmin(wsp)), c(ymax(wsp), xmin(wsp)), fun = distHaversine) # coordinat system
ydist=distm (c(ymin(wsp), xmin(wsp)), c(ymin(wsp), xmax(wsp)), fun = distHaversine)
stepx=2 #step in meters axis X
stepy=2 #step in meters axis Y
ncolgr=as.vector(xdist%/%stepx)
nrowgr=as.vector(ydist%/%stepy)

xy.sp=SpatialPoints(dataRMO[1:2])
xy.spdf<-SpatialPointsDataFrame(xy.sp,log(dataRMO[3])) # interpolation of logarithm

xg=seq(xmin(wsp),xmax(wsp), length.out=ncolgr) #step in meters
yg=seq(ymin(wsp),ymax(wsp), length.out=nrowgr)

pts=expand.grid(xg,yg) # grid for results of interpolation 
grd.pts=SpatialPixels(SpatialPoints(pts))
proj4string(xy.spdf)=CRS("+init=epsg:4326")
proj4string(grd.pts)=CRS("+init=epsg:4326")

x <- krige(z~1, xy.spdf, grd.pts) #interpolation
OutGrid=cbind(x$Var1,x$Var2,exp(x$var1.pred)) # exponentiation
df1=as.data.frame(OutGrid)
write.table(df1,file=paste("GRidded_",inputdatafilename,sep=""),sep=" ")
