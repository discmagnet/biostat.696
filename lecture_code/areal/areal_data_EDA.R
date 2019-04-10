

###  This file shows how to do EDA for areal data
###  The packages that we will use for this purpose are
###  maps, maptools, spdep, classInt, RColorBrewer and SpatialEpi

library(maps)
library(maptools)
library(spdep)
library(classInt)
library(RColorBrewer)
library(SpatialEpi)


### The maps package is basic a database of shapefiles for the world, the US and some foreign countries
### To see some of the database that are contained in maps, type ?map

? map
## Here we are making a map of the world with country boundaries
map("world",interior=TRUE,plot=TRUE)
## Here we are saving the database world in an R object called wrld
wrld <- map("world",fill=TRUE,plot=FALSE)
names(wrld)
## This is the name of all the polygons that are contained in the world database
wrld.names <- wrld$names
length(wrld.names)

###  Similarly, we can look at the US

## this draws the map of the US states with boundaries among the states
map("state",interior=TRUE)
## this draws the map of the US states without boundaries among the states
map("state",interior=FALSE)

## this is to save the object US that is stored in the database maps
## and pass it to the sp and spdep packages for computation purposes
US <- map("state",fill=TRUE,plot=FALSE)
US
names(US)

range(US$x,na.rm=TRUE)
range(US$y,na.rm=TRUE)

US.names <- US$names
length(US$names)

## this is to get IDs that are useful later on
US.IDs <- sapply(strsplit(US.names,":"),function(x) x[1])

## Here we transform the US dataset in the maps database into a SpatialPolygon object
US_poly_sp <- map2SpatialPolygons(US,IDs=US.IDs,proj4string=CRS("+proj=longlat + datum=wgs84"))
plot(US_poly_sp,col="white",axes=TRUE)

## To know the names of the states as they appear in the US_poly_sp object, we do
row.names(US_poly_sp)
length(row.names(US_poly_sp))

## If for example, we dont want the district of columbia in our dataset,
## we can modify the spatial polygon object and remove from it
## the row that refers to district of columbia
index.columbia <- 8
US_no_columbia <- US_poly_sp[seq(1,49)[-8]]

row.names(US_no_columbia)

## This just draws the border of the US without the district of Columbia
plot(US_no_columbia,border="black",axes=TRUE)

## We can access the information regarding each state in the US_no_columbia 
## by simply doing US_no_columbia[#i] where i is a number between 1 and 48
US_poly_sp[1]
US_poly_sp[48]

## We can color Alabama blue, Arizona red, Georgia green and Florida orange 
## by simply doing the following
plot(US_no_columbia[1],col="red",add=T)
plot(US_no_columbia[2],col="blue",add=T)
plot(US_no_columbia[8],col="orange",add=T)
plot(US_no_columbia[9],col="green",add=T)


#####

## Reading the geographical information from a shape file
## Here we are reading a dataset that contains information on the percentage of sample with blood 
## group A for given towns/area in Eire
## Other information in the dataset include: the percentage of people in the town/area
## that is beyond pale and people whose skin complexion if pale; the number of blood type samples;
## the arterial road network accessability in 1961, the percentage in value terms of 
## gross arterial output of each county consumed by itself, the population in 1961 as a
## percentage of the population in 1926, the value of retail sales in pound,
## the total personal income in pounds, and the county names
eire <- readShapePoly(system.file("etc/shapes/eire.shp", package="spdep")[1],
ID="names", proj4string=CRS("+proj=utm +zone=30 +units=km"))

names(eire)

row.names(eire)

A.blood <- eire$A
towns <- eire$towns
pale <- eire$pale
size <- eire$size
income <- eire$INCOME
sale <- eire$RETSALE
eire.county <- eire$names


##############  Clorophleth maps

### Here we make a plot of the SIDS rate in North Carolina during the period July 1, 1974 to June 30, 1978
nc.sids <- read.csv("SIDS1974_NC.csv",sep=",")
dim(nc.sids)
[1] 100   4

names(nc.sids)
[1] "County"           "No_births"        "No_deaths"        "Non.white_births"

county <- nc.sids$County
no.births <- nc.sids$No_births
no.deaths <- nc.sids$No_deaths
nw.births <- nc.sids$Non.white_births


## Here we are getting the information from the maps package in R about the counties in NC
nc.counties <-map("county","north carolina",fill=T,plot=F)
## This is to get the IDs of the counties in NC that are useful to make the spatial polygon
nc.IDs <- sapply(strsplit(nc.counties$names,","),function(x) x[2])
nc.poly <- map2SpatialPolygons(nc.counties,IDs=nc.IDs,proj4string=CRS("+proj=longlat + datum=wgs84"))

## This command tells R to derive the neighborhood structure for the counties in North Carolina
nc.nb <- poly2nb(nc.poly)
nc.nb

## Here we get the adjacency weights, together with the information on the number of weights for each county
us.weights <- nb2WB(us.nb)
adj.us <- us.weights$adj
weights.us <- us.weights$weights
num.us <- us.weights$num
us.listw <- nb2listw(us.nb)



plotvar <- A.blood
nclr <- 9

plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))
class

### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)


plot(eire,border="black",axes=TRUE,xlim=c(-200,400))
title(xlab="Eastings (km)",ylab="Northings (km)",main="Percentage of blood A samples")
plot(eire,col=colcode,add=T)

leg.txt<-c("[23.9%,25.2%)","[25.2%,26.6%)","[26.6%,27.9%)","[27.9%,29.2%)",
"[29.2%,30.6%)","[30.6%,31.9%)","[31.9%,33.2%)","[33.2%,34.5%)","[34.5%,35.9%]")
legend(-200,6100,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

#### Using a different definition of class intervals
plotvar <- A.blood
nclr <- 4

plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="quantile")
class

colcode <- findColours(class,plotclr)


plot(eire,border="black",axes=TRUE,xlim=c(-200,400))
title(xlab="Eastings (km)",ylab="Northings (km)",main="Percentage of blood A samples")
plot(eire,col=colcode,add=T)

leg.txt<-c("[23.9%,27.9%)","[27.9%,29.3%)","[29.3%,30.9%)","[30.9%,35.9%]")
legend(-200,6100,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

#### Using a different definition of class intervals
plotvar <- A.blood
nclr <- 9

plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="equal")
class

colcode <- findColours(class,plotclr)

plot(eire,border="black",axes=TRUE,xlim=c(-200,400))
title(xlab="Eastings (km)",ylab="Northings (km)",main="Percentage of blood A samples")
plot(eire,col=colcode,add=T)

leg.txt<-c("[23.9%,25.2%)","[25.2%,26.6%)","[26.6%,27.9%)","[27.9%,29.2%)",
"[29.2%,30.5%)","[30.5%,31.9%)","[31.9%,33.2%)","[33.2%,34.5%)","[34.5%,35.9%]")
legend(-200,6100,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")




### Here we make a plot of the observed income
plotvar <- income
nclr <- 9

plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))
class

### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)


plot(eire,border="black",axes=TRUE,xlim=c(-210,400))
title(xlab="Eastings (km)",ylab="Northings (km)",main="Income")
plot(eire,col=colcode,add=T)

leg.txt<-c("[5,297;23,000.8%)","[23,000.8;40,704.6)","[40,704.6;58,408.3)","[58,408.3;76,112.1)",
"[76,112.1;93,815.9)","[93,815.9;111,519.7)","[111,519.7;129,223.4)","[129,223.4;146,927.2)","[146,927.2;164,631.0]")
legend(-210,6100,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")



##############  Building the spatial proximity matrix

##### Neighbors based on the definition of contiguity

## Here we can use the spdep package to get the list of neighbors 
## from the spatial polygon automatically
## The default definition of neighbors is the "Queen" definition
## If the option queen is set equal to FALSE, then we are using the "Rook"
## definition of neighbors
eire.nb <- poly2nb(eire)
eire.nb

### This makes a map of Ireland with a link between counties that are adjacent according 
### to the definition of neighbors that we have defined
### This first command takes the coordinate of the centroid for each of the 26 counties in Ireland
eire.coord <- coordinates(eire)
plot(eire,border="black",axes=TRUE)
### Here we add the link between states that are adjacent
plot(eire.nb,eire.coord,col="red",pch=20,add=TRUE)


##### K nearest neighbors

## The function knearneighb returns a matrix with the indices of the k nearest neighbors to
## each areal unit
eire.knn <- knearneigh(eire.coord,k=2)
## If we apply to it the function knn we get an object of the type nb to which we can apply the function
## nb2WB
eire.knn.nb <- knn2nb(eire.knn,row.names=eire.county)
eire.knn.nb

eire.knn.weights <- nb2WB(eire.knn.nb)
eire.knn.weights
adj.eire.knn <- eire.knn.weights$adj
weights.eire.knn <- eire.knn.weights$weights
num.eire.knn <- eire.knn.weights$num


## Here we make a map of Ireland with a link between counties that are adjacent according to this 
## definition of neighbors
eire.coord <- coordinates(eire)
plot(eire,border="black",axes=TRUE)
### Here we add the link between states that are adjacent
plot(eire.knn.nb,eire.coord,col="red",pch=20,add=TRUE)


##### Neighbord based on thresholding the distance

## Here we determine the neighbors of each of the 26 counties in Ireland using a definition
## based on thresholding the distance among units

## Here we compute all the distances among the representative for each county
eire.dsts <- unlist(nbdists(eire.nb,eire.coord,longlat=FALSE))

## Then we look at the distribution of these distances, and we decide that we
## define neighbors areal units for which the distance is greater than 0 and at more 
## 75% of the maximum distance
## The function that automatically generates the weights according to this definition
## is the function dnearneigh
summary(eire.dsts)
eire.n.tresh1 <- dnearneigh(eire.coord,d1=0,d2=0.75*max(eire.dsts),longlat=FALSE)

## Then, we can use the function nb2WB as usual to obtain an object that is more easily
## interpretable and that can be mapped
eire.n.tresh1.weights <- nb2WB(eire.n.tresh1)
eire.n.tresh1.weights

adj.eire.n.tresh1 <- eire.n.tresh1.weights$adj
weights.eire.n.tresh1 <- eire.n.tresh1.weights$weights
num.eire.n.tresh1 <- eire.n.tresh1.weights$num


## Here we make a map of Ireland with a link between counties that are adjacent according to this 
## definition of neighbors
eire.coord <- coordinates(eire)
plot(eire,border="black",axes=TRUE)
### Here we add the link between states that are adjacent
plot(eire.n.tresh1,eire.coord,col="red",pch=20,add=TRUE)


##############  Computing Morans I

## To compute Moran I, we need a set of weights
## In spdep the weights that need to be passed to the moran.test function needs to be in a list form
## So, we apply the nb2listw function on the nb object to get a set of weigths in a list form
eire.list.w <- nb2listw(eire.nb)

## The function moran.test computes Moran I for areal data and returns the observed
## value of I, the expected value under the hypothesis of iid data and 
## also the asymptotic variance (we saw the formula in class)
## If we specify randomisation=FALSE, the function moran.test also
## runs a hypothesis test on I using the asymptotic distribution of I
moran.eire <- moran.test(A.blood,listw=eire.list.w,randomisation=FALSE)
moran.eire

## Here we are computing Moran I on the residuals of a linear regression of the percentage of A blood samples
## and the indiator for pale or not pale
lm.morantest(lm(A.blood~pale),listw=eire.list.w)

## The function moran.mc runs a hypothesis test on Moran if
## using a Monte Carlo approach
## For this reason we need to specify the number of samples to simulate from all the possible N! permutations
## of the data
bperm.moran.eire <- moran.mc(A.blood,listw=eire.list.w,nsim=9999)
bperm.moran.eire

## Here we are running a hypothesis test on Moran I for the residuals of a linear regression of the percentage of A blood samples
## and the indiator for pale or not pale using a Monte Carlo approach
res.moran.eire <- as.numeric(lm(A.blood~pale)$residuals)
bperm.res.moran.eire <- moran.mc(res.moran.eire,listw=eire.list.w,nsim=9999)
bperm.res.moran.eire


##############  Computing Geary C

## The syntax for this funtion is exactly the same as for Moran I
## So, we apply the nb2listw function on the nb object to get a set of weigths in a list form

eire.list.w <- nb2listw(eire.nb)
geary.eire <- geary.test(A.blood,listw=eire.list.w,randomisation=FALSE)
geary.eire


bperm.geary.eire <- geary.mc(A.blood,listw=eire.list.w,nsim=9999)
bperm.geary.eire

res.geary.eire <- as.numeric(lm(A.blood~pale)$residuals)
bperm.res.geary.eire <- geary.mc(res.geary.eire,listw=eire.list.w,nsim=9999)
bperm.res.geary.eire


############# Correlogram

## The function sp.correlogram computes the correlogram for the spatial variable.
## We can specify the maximum order of neighbors that we want to consider, what index
## of spatial autocorrelation to use (choices are: Moran I and Geary C)
## and we need to specify a policy for how to handle cases in which 
## an area has no neighbors. The default is zero.policy=FALSE
## in which case the function gives an error if there are areal units with
## no neighbors. Specifying zero.policy=T avoids this problem

cor.mi.5 <- sp.correlogram(neighbours=mi.nb,var=unempt,
order=5,method="I",style="C",zero.policy=T)
print.spcor(cor.mi.5)

## This is to make a plot of the correlogram
plot.spcor(cor.mi.5,main="Correlogram for unemployment rate")


############# Algorithm for spatial filtering

## Here we show how to
cor.mi.5 <- sp.correlogram(neighbours=mi.nb,var=unempt,
order=5,method="I",style="C",zero.policy=T)
print.spcor(cor.mi.5)

## This is to make a plot of the correlogram
plot.spcor(cor.mi.5,main="Correlogram for unemployment rate")


