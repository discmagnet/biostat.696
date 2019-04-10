setwd("~/WORKING_DIRECTORIES/biostat.696/lecture_code")

####################################################################################
#####
#####  This code illustrates how to make spatial maps of spatial data in R
#####
####################################################################################

##############################
# Note. To install a packge, you 
# first need to download the package
# Once it is downloaded and installed 
# (for help on installing packages in R, 
# type ?INSTALL at the prompt)
# you can upload it using the command 
# library(name of the package). 
# An example is given below.
#
# To unload a package use the command
# detach("package:name of the package")
##############################

# Point-referenced data

# We will use the ozone dataset that you can find on Ctools in the Data folder.
# The file is named: "Ozone_Aug1_2001.txt"

# Here we read the data in
ozone <- read.table("Ozone_Aug1_2001.txt",header=TRUE)
# This is to obtain more information on the possible options that can be specified with the command "read.table"
?read.table

# This command is to see how many lines and variables there are in the file
dim(ozone)

# The file should contain 1128 ozone measurements and 4 variables: 
# the site latitude, the site longitude, the date of observation and the observed ozone concentration at the site.
# Here we check if we have read in the dataset correctly
ozone[1:3,]

lon <- ozone$longitude
lat <- ozone$latitude
o3 <- ozone$ozone


plot(lon,lat,xlab="Longitude",ylab="Latitude",type="p",pch=20,col="black")

library(fields)
plot(lon,lat,xlab="Longitude",ylab="Latitude",type="p",pch=20,col="black")
US(xlim=range(lon),ylim=range(lat),add=T,lwd=1,col="gray")

# If we want to change the range for the US map, make sure to change the limits also in the plot command line.
plot(lon,lat,xlab="Longitude",ylab="Latitude",type="p",pch=20,col="black",xlim=c(-130,-65),ylim=c(22,52))
US(xlim=c(-130,-65),ylim=c(22,52),add=T,lwd=1,col="gray")


# The commands above just map the locations of the sites with observations but not the observed values
# The RgoogleMaps allows to point the observations locations on a Google Map
library(RColorBrewer)
library(classInt)
library(RgoogleMaps)

# Here we show how to overlay the observation locations on a Google Map
MyMap <- GetMap.bbox(lonR=range(lon),latR=range(lat),size=c(640,640),maptype="hybrid")
PlotOnStaticMap(MyMap)

# See book on p.18 for the explanation of what this command is doing
convert_points <- LatLon2XY.centered(MyMap,lat,lon)
points(convert_points$newX,convert_points$newY,col="red",pch=19)


# Here we define how many different colors we want to use in the map and the variable to plot
# This is variable to plot
plotvar <- o3
# This is the number of colors that we want to use (Note that the number of colors that can be used depend on the 
# palette of colors chosen. See ?brewer.pal for more information)
nclr <- 9
plotclr <- brewer.pal(nclr,"YlOrRd")
# This is a different way to create a color palette with 10 different colors
col.scheme <- colorRampPalette(c("blue","cyan","yellow","red"))
plotclr.2 <- col.scheme(nclr)

# This command divides the variable to plot in nclr classes with fixed breaks. The break points are specified in
# the fixedBreaks. Other methods can be chose. See ?classIntervals
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))

#class is a list with two variables: "var" and "breaks"
names(class)
# These are the breakpoints of the class variable
class$brks
# This command assigns a color to each class. If we want to use the second color scheme, we change plotclr to plotclr.2
colcode <- findColours(class,plotclr)


plot(lon,lat,xlab="Longitude",ylab="Latitude",type="p",pch=20,col=colcode,main="Ozone concentration on 08/31/2001")
US(xlim=range(lon),ylim=range(lat),add=T,lwd=1,col="gray")


# To save the plot as a postscript file
postscript("ozone_concentration.ps",height=9,width=9)
plot(lon,lat,xlab="Longitude",ylab="Latitude",type="p",pch=20,col=colcode,main="Ozone concentration on 08/31/2001")
US(xlim=range(lon),ylim=range(lat),add=T,lwd=1,col="gray")
dev.off()


# To add a legend, we need to use a different function: the function image.plot
# The data has to be given in a specific format: see the help files for 
# image.plot or better image (in the graphics package)
?image.plot
?image

# The function image requires that the x and y coordinates be 
# grids in ascending order of the same length and that the z vector be a 
# matrix of the same dimension as the x and y coordinate vectors

length.x <- 25
length.y <- 25

x.grid <- seq(min(lon),max(lon),length=length.x)
y.grid <- seq(min(lat),max(lat),length=length.y)
z.grid <- matrix(NA,nrow=length.x,ncol=length.y)
z.lim <- range(plotvar)

# This first command does not really plot anything except for the legend
image.plot(x.grid,y.grid,z.grid,xlab="Latitude",ylab="Longitude",zlim=z.lim,
col=plotclr,breaks=class$brks,legend.lab="ppb",main="Ozone concentration on 08/31/2001")
US(xlim=range(lon),ylim=range(lat),add=T,lwd=1,col="gray")
points(lon,lat,col=colcode,pch=20)



# This command is to make a map of a surface, obtained by interpolating the
# observed values onto a grid
library(akima)
surf.ozone <- interp(lon,lat,o3)

# We can make a plot of this surface directly
image.plot(surf.ozone$x,surf.ozone$y,surf.ozone$z,col=plotclr,breaks=seq(min(surf.ozone$z,na.rm=T),max(surf.ozone$z,na.rm=T),length=10),
zlim=range(surf.ozone$z,na.rm=T),xlab="Longitude",ylab="Latitude",legend.lab="ppb",main="Interpolated ozone concentration on 08/31/2001")
US(xlim=range(lon),ylim=range(lat),add=T,lwd=1,col="black")


# This is using a different color scheme
image.plot(surf.ozone$x,surf.ozone$y,surf.ozone$z,col=heat.colors(100)[100:1],
zlim=range(surf.ozone$z,na.rm=T),xlab="Longitude",ylab="Latitude",legend.lab="ppb",
main="Interpolated ozone concentration on 08/31/2001")
US(xlim=range(lon),ylim=range(lat),add=T,lwd=2,col="black")
# To add contour lines
contour(surf.ozone$x,surf.ozone$y,surf.ozone$z,nlevel=10,add=T,col="black")

# This is to make a 3-D plot of the interpolated surface
drape.plot(surf.ozone$x,surf.ozone$y,surf.ozone$z,xlab="Longitude",ylab="Latitude",zlab="ppb",
col=heat.colors(100)[100:1],theta=15,phi=20,border=FALSE,main="Interpolated ozone concentration on 08/31/2001")



##############################

# This is to show how to compute great-circle distances in R using the package fields
rm(list=ls())
library(fields)

ozone <- read.table("Ozone_Aug1_2001.txt",header=TRUE)
lon <- ozone$longitude
lat <- ozone$latitude
o3 <- ozone$ozone


# Here we create a matrix with two columns: the first column contains the longitude of the monitoring sites
# the second column contains the latitude of the monitoring sites
matrix.sites <- matrix(cbind(lon,lat),nrow=length(lon),ncol=2)
matrix.sites[1:3,]

# The function in the fields package that measures the great circle distance among points on the Earth is rdist.earth
# The output is a square symmetric matrix whose (i,j) element is the distance between site_i and site_j
Distance.matrix <- rdist.earth(matrix.sites,matrix.sites,miles=FALSE)

Distance.matrix[1,1]
Distance.matrix[1,2]
Distance.matrix[1,3]



## This it to compute the distance between points using Euclidean distance
## and the 3-dimensional coordinates 

## Here we are deriving the 3-D coordinates of each point
R <- 6371
x.proj <- R*cos((lat*pi)/180)*cos((lon*pi)/180)
y.proj <- R*cos((lat*pi)/180)*sin((lon*pi)/180)
z.proj <- R*sin((lat*pi)/180)
matrix.sites.proj <- matrix(cbind(x.proj,y.proj,z.proj),nrow=length(x.proj),ncol=3)

## Here we are computing the euclidean distance between points
Distance.matrix.eucl <- rdist(matrix.sites.proj,matrix.sites.proj)
Distance.matrix.eucl[1,1]
Distance.matrix.eucl[1,2]
Distance.matrix.eucl[1,3]


## This is to show how to convert longitude and latitude into UTM coordinates. To do this, we can use the rgdal package.
## Note that in order to be able to convert coordinates into UTM coordinate, it is necessary to specify the UTM zone explicitly

# We are using data collected in Colorado
coloradoST <- read.table("ColoradoS-T.dat",header=TRUE)
library(sp)
library(rgdal)

SP_longlat <- SpatialPoints(coords=cbind(coloradoST$Lon,coloradoST$Lat),proj4string=CRS("+proj=longlat +ellps=WGS84"))
# This transforms longitude and latitude in UTM coordinates
SP_utm <- spTransform(SP_longlat,CRS("+proj=utm +zone=13 +datum=WGS84"))




