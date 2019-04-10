setwd("~/WORKING_DIRECTORIES/biostat.696/lecture_code")

####################################################################################
#####
#####  This code illustrates how to run EDA for spatial data in R
#####
####################################################################################

library(MASS)
library(fields)
library(akima)
library(RColorBrewer)
library(classInt)
library(gstat)
library(sp)


###### Mean/median versus geographical coordinate

# We will use the soil dataset that you can find on Ctools in the Data folder
# The file is names: "Soil_data.txt"

soil <- read.table("Soil_data.txt",sep=" ",header=TRUE)

# These are the geographical coordinates
x.soil <- soil$x
y.soil <- soil$y
# This our spatial process
water.soil <- soil$pH_water


plotvar <- water.soil
nclr <- 9
plotclr <- brewer.pal(nclr,"Blues")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))

names(class)
class$brks
colcode <- findColours(class,plotclr)

length.x <- 25
length.y <- 25

x.grid <- seq(min(x.soil),max(x.soil),length=length.x)
y.grid <- seq(min(y.soil),max(y.soil),length=length.y)
z.grid <- matrix(NA,nrow=length.x,ncol=length.y)
z.lim <- range(plotvar)

# This first command does not really plot anything except for the legend
image.plot(x.grid,y.grid,z.grid,xlab="X",ylab="Y",zlim=z.lim,
col=plotclr,breaks=class$brks,legend.lab="pH",main="Water pH")
points(x.soil,y.soil,col=colcode,pch=20)


# Here we take the unique values for the x (resp. y) coordinate and we sort them in ascending order
un.x.soil <- sort(unique(x.soil))
un.y.soil <- sort(unique(y.soil))

# Here we create the vectors that will contain the mean and median values of our
# spatial variable, the Ph content of water for each unique value of the x and the y coordinate
# that is for each row and column
mean.water.x <- rep(0,length(un.x.soil))
mean.water.y <- rep(0,length(un.y.soil))
median.water.x <- rep(0,length(un.x.soil))
median.water.y <- rep(0,length(un.y.soil))

# This is the loop to fill each entry of these vectors
for(i in 1:length(un.x.soil)){
	mean.water.x[i] <- mean(water.soil[x.soil==un.x.soil[i]])	
	median.water.x[i] <- quantile(water.soil[x.soil==un.x.soil[i]],0.5)	
}

for(i in 1:length(un.y.soil)){
	mean.water.y[i] <- mean(water.soil[y.soil==un.y.soil[i]])	
	median.water.y[i] <- quantile(water.soil[y.soil==un.y.soil[i]],0.5,na.rm=T)	
}

# Here we make plots of the mean/median of the spatial variable versus the unique values of the x (resp. y) coordinate
plot(un.x.soil,mean.water.x,type="p",col="black",pch=20,
ylim=range(c(mean.water.x,median.water.x)),main="Plot of the mean/median of Y(s_i) versus the x coordinate",
xlab="X",ylab="pH of Water")
points(un.x.soil,median.water.x,pch=20,col="red")
legend("topright",pch=rep(20,2),col=c("black","red"),c("Mean","Median"))

plot(un.y.soil,mean.water.y,type="p",col="black",pch=20,
ylim=range(c(mean.water.y,median.water.y)),main="Plot of the mean/median of Y(s_i) versus the y coordinate",
xlab="Y",ylab="pH of Water")
points(un.y.soil,median.water.y,pch=20,col="red")
legend("topright",pch=rep(20,2),col=c("black","red"),c("Mean","Median"))


#####  Boxplots versus index or row column

# Here we make boxplots of the spatial variable versus the x (resp.y) coordinate
# First we need to define a factor variable for the x and y coordinate
grp.x <- as.factor(x.soil)
grp.y <- as.factor(y.soil)

# This makes the boxplot of water ph grouped with respect to the x coordinate
boxplot(water.soil~grp.x,col="grey",xlab="X",ylab="pH of Water")

# Here is with respect to the y coordinate
boxplot(water.soil~grp.y,col="grey",xlab="Y",ylab="pH of Water")


#### Estimating spatial trend

x.soil.sq <- x.soil*x.soil
y.soil.sq <- y.soil*y.soil
x.y.soil <- x.soil*y.soil
x.soil.cub <- x.soil*x.soil*x.soil
x.soil.sq.y.soil <- x.soil*x.soil*y.soil
x.soil.y.soil.sq <- x.soil*y.soil*y.soil
y.soil.cub <- y.soil*y.soil*y.soil

# Here we fit a trend surface model: a full third degree polynomial in the geographical coordinates 
mod.water <- lm(water.soil~x.soil+y.soil+x.soil.sq+y.soil.sq+(x.y.soil)+x.soil.cub+x.soil.sq.y.soil+x.soil.y.soil.sq+y.soil.cub)
coeff.water <- summary(mod.water)$coeff[1:10]

# The estimated spatial trend then is given by the fitted values and we can obtain a surface by bilinarly interpolating
# using the R command interp
trend.water <- mod.water$fitted.values
surf.trend.water <- interp(x.soil,y.soil,trend.water)

# We want to use the same classes and breaks as the original data
image.plot(surf.trend.water$x,surf.trend.water$y,surf.trend.water$z,col=plotclr,
breaks=class$brks,zlim=range(water.soil,na.rm=T),xlab="X",ylab="Y",main="Estimated spatial trend")


################################################
###### Variogram 
################################################

# To illustrate the example with the empirical variogram we use a different dataset
# We will use the meuse dataset that is available on the gstat package and also saved on the Ctools website as the "meuse.txt" file

library(fields)
library(gstat)

data(meuse)
x.zinc <- meuse$x
y.zinc <- meuse$y
zinc <- meuse$zinc

# We work with the log zinc concentration, because the distribution of zinc is skewed and to reduce the variability
l.zinc <- log(zinc)

# Here we estimate the trend using a quadratic trend surface model
x.zinc.sq <- x.zinc*x.zinc
y.zinc.sq <- y.zinc*y.zinc
x.y.zinc <- x.zinc*y.zinc

mod.lzinc.sq <- lm(l.zinc~x.zinc+y.zinc+x.zinc.sq+y.zinc.sq+x.y.zinc)
coeff.lzinc.sq <- summary(mod.lzinc.sq)$coeff[1:6]
trend.lzinc <- mod.lzinc.sq$fitted.values 
res.lzinc.sq <- l.zinc-trend.lzinc


surf.est.surface.lzinc.sq <- interp(x.zinc,y.zinc,trend.lzinc)

plotvar <- surf.est.surface.lzinc.sq$z[!is.na(surf.est.surface.lzinc.sq$z)]
nclr <- 9
plotclr <- brewer.pal(nclr,"YlOrRd")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))

names(class)
class$brks
colcode <- findColours(class,plotclr)



image.plot(surf.est.surface.lzinc.sq$x,surf.est.surface.lzinc.sq$y,
surf.est.surface.lzinc.sq$z,col=plotclr,breaks=seq(4,8,length=10),
zlim=c(4,8),xlab="X",ylab="Y",main="Estimated spatial trend")


## Variogram cloud plot

# Here we set up the distance matrix with all the distances among sites
n.zinc <- length(x.zinc)
Dist.mat <- matrix(0,nrow=n.zinc,ncol=n.zinc)

coord.zinc <- matrix(cbind(x.zinc,y.zinc),nrow=n.zinc,ncol=2)
# The rdist function computes the euclidean distance
Dist.mat <- rdist(coord.zinc,coord.zinc)

# Here we set up the matrix with the squared differences
Diff.mat <- matrix(0,n.zinc,n.zinc)

for(i in 1:n.zinc){
	for(j in 1:n.zinc){
		Diff.mat[i,j] <- (1/2)*(res.lzinc.sq[i]-res.lzinc.sq[j])^2	
	}
}


# This is a plot of the variogram cloud
plot(as.numeric(Dist.mat),as.numeric(Diff.mat),type="p",pch=20,col="black",xlab="Distance",ylab="Squared differences")

# To compute the empirical variogram, we use the gstat package
# First we need to create a dataframe with the data
lzinc.data <- data.frame(cbind(x.zinc,y.zinc,res.lzinc.sq))

emp.variog.lzinc <- variogram(res.lzinc.sq~1,locations=~x.zinc+y.zinc,lzinc.data)
emp.variog.lzinc
# This is to plot the empirical variogram
plot(emp.variog.lzinc,col="black",type="p",pch=20)



## This is to compute the empirical directional variograms
dir.variog.lzinc <- variogram(res.lzinc.sq~1,locations=~x.zinc+y.zinc,lzinc.data,alpha=c(0,45,90,135))
# This command plots each empirical directional variogram in a separate frame
plot(dir.variog.lzinc)

# To plot each of them in the same frame and have them superimposed, we need to first identify the values corresponding to each angle
dir.variog.lzinc.E <- dir.variog.lzinc[dir.variog.lzinc$dir.hor==0,]
dir.variog.lzinc.NE <- dir.variog.lzinc[dir.variog.lzinc$dir.hor==45,]
dir.variog.lzinc.N <- dir.variog.lzinc[dir.variog.lzinc$dir.hor==90,]
dir.variog.lzinc.NW <- dir.variog.lzinc[dir.variog.lzinc$dir.hor==135,]

plot(dir.variog.lzinc.E$dist,dir.variog.lzinc.E$gamma,col="blue",type="o",pch=20,ylim=c(0,1),xlab="Distance",ylab="Semivariance")
lines(dir.variog.lzinc.NE$dist,dir.variog.lzinc.NE$gamma,col="orange",type="o",pch=20)
lines(dir.variog.lzinc.N$dist,dir.variog.lzinc.N$gamma,col="green",type="o",pch=20)
lines(dir.variog.lzinc.NW$dist,dir.variog.lzinc.NW$gamma,col="magenta",type="o",pch=20)
legend(50,1,col=c("blue","orange","green","magenta"),lwd=rep(1,5),lty=rep(1,5),
c("E-W","NE-SW","N-S","NW-SE"))

