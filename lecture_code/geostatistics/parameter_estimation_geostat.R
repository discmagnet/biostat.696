setwd("~/WORKING_DIRECTORIES/biostat.696/lecture_code")

####################################################################################
#####
#####  This code illustrates how to estimate parameters for spatial processes in R
#####
####################################################################################

library(MASS)
library(fields)
library(akima)
library(RColorBrewer)
library(classInt)
library(gstat)


# We will use the soil dataset that you can find on Ctools in the Data folder
# The file is names: "Soil_data.txt"

soil <- read.table("Soil_data.txt",sep=" ",header=TRUE)

# These are the geographical coordinates
x.soil <- soil$x
y.soil <- soil$y
# This our spatial process
water.soil <- soil$pH_water


## This is the part where we estimate the spatial trend
## We are using a surface trend model with a 3rd degree polynomial
x.soil.sq <- x.soil*x.soil
y.soil.sq <- y.soil*y.soil
x.y.soil <- x.soil*y.soil
x.soil.cub <- x.soil*x.soil*x.soil
x.soil.sq.y.soil <- x.soil*x.soil*y.soil
x.soil.y.soil.sq <- x.soil*y.soil*y.soil
y.soil.cub <- y.soil*y.soil*y.soil

mod.water <- lm(water.soil~x.soil+y.soil+x.soil.sq+y.soil.sq+(x.y.soil)+x.soil.cub+x.soil.sq.y.soil+x.soil.y.soil.sq+y.soil.cub)
coeff.water <- summary(mod.water)$coeff[1:10]

# These are the OLS estimates of the coefficients
coeff.water

# These are the fitted values 
trend.water <- mod.water$fitted.values

# These are the residuals
res.water <- water.soil-trend.water

# Here we are computing the empirical variogram using the variogram function in gstat
# In order to do that, we need to put the data in the right format
# So, we define a data frame, and we specify which ones are the coordinates in the data.frame
water.data <- data.frame(cbind(x.soil,y.soil,res.water))
emp.variog.water <- variogram(res.water~1,locations=~x.soil+y.soil,water.data)

plot(emp.variog.water)


############################################################
#
#     WLS estimation of the variogram parameters
#
############################################################

## To fit a parametric variogram model via WLS we use the fit.variogram function in the gstat
## library

## To see what are the parametric variograms implemented in the gstat library, look at the vgm function
?vgm

## We will try: Exponential, Spherical, Gaussian and Matern
## Inspecting the empirical variogram plot, it seems as if there is a nugget effect, equal to 0.01
## the sill seems to be 0.02 and the partial sill is therefore 0.01
## the range seems to be 10

# Here we try the exponential semi-variogram
# Note that we specify fit.method=2 to do WLS. The default for fit.method is to use weights equal to N(h_m)/((h_m)^2)
exp.variog.water <- fit.variogram(emp.variog.water,vgm(psill=0.01,"Exp",10,0.01),fit.method=2)
exp.variog.water
model      psill    range
1   Nug 0.00000000 0.000000
2   Exp 0.02175365 5.730311

# This is to make plots of the empirical variogram with overlaid the exponential variogram
print(plot(emp.variog.water,exp.variog.water))

# This is to do it by hand:
exp.variog <- function(sigma2,phi,tau2,dist){
	n <- length(dist)
	exp.variog.vec <- rep(0,n)
	for(i in 1:n){
		exp.variog.vec[i] <- tau2+(sigma2*(1-exp(-dist[i]/phi)))
	}
	return(exp.variog.vec)
}


sigma2.exp <- exp.variog.water$psill[2]
phi.exp <- exp.variog.water$range[2]
tau2.exp <- exp.variog.water$psill[1]
dist.vec <- emp.variog.water$dist


exp.variog.water.2 <- exp.variog(sigma2.exp,phi.exp,tau2.exp,dist.vec)

plot(dist.vec,emp.variog.water$gamma,type="p",pch=20,col="black",xlab="Distance",ylab="Semi-variance",main="Exponential semi-variogram",ylim=c(0,0.025))
lines(dist.vec,exp.variog.water.2,col="red")


# Here we try the spherical semi-variogram
sph.variog.water <- fit.variogram(emp.variog.water,vgm(psill=0.01,"Sph",range=10,nugget=0.01),fit.method=2)
sph.variog.water
model       psill    range
1   Nug 0.003396606  0.00000
2   Sph 0.017935018 15.49115

# This is to make plots of the empirical variogram with overlaid the spherical variogram
print(plot(emp.variog.water,sph.variog.water))

# This is to do it by hand:
sph.variog <- function(sigma2,phi,tau2,dist){
	n <- length(dist)
	sph.variog.vec <- rep(0,n)
	for(i in 1:n){
		if(dist[i] < phi){
			sph.variog.vec[i] <- tau2+(sigma2*(((3*dist[i])/(2*phi))-((dist[i]^3)/(2*(phi^3)))))
		}
		if(dist[i] >= phi){
			sph.variog.vec[i] <- tau2+sigma2	
		}
	}
	return(sph.variog.vec)
}


sigma2.sph <- sph.variog.water$psill[2]
phi.sph <- sph.variog.water$range[2]
tau2.sph <- sph.variog.water$psill[1]
dist.vec <- emp.variog.water$dist

sph.variog.water.2 <- sph.variog(sigma2.sph,phi.sph,tau2.sph,dist.vec)

plot(dist.vec,emp.variog.water$gamma,type="p",pch=20,col="black",xlab="Distance",ylab="Semi-variance",main="Spherical semi-variogram",ylim=c(0,0.025))
lines(dist.vec,sph.variog.water.2,col="red")


# Here we try the Gaussian semi-variogram
gau.variog.water <- fit.variogram(emp.variog.water,vgm(psill=0.01,"Gau",range=10,nugget=0.01),fit.method=2)
gau.variog.water
model       psill    range
1   Nug 0.007250089 0.000000
2   Gau 0.014050270 8.357302


sigma2.gau <- gau.variog.water$psill[2]
phi.gau <- gau.variog.water$range[2]
tau2.gau <- gau.variog.water$psill[1]


# Here we try the Matern semi-variogram
mat.variog.water <- fit.variogram(emp.variog.water,vgm(psill=0.01,"Mat",range=10,nugget=0.01,kappa=1),fit.method=2)
mat.variog.water
model      psill  range kappa
1   Nug 0.00000000 0.0000     0
2   Mat 0.02167852 3.6725     1

sigma2.mat <- mat.variog.water$psill[2]
phi.mat <- mat.variog.water$range[2]
tau2.mat <- mat.variog.water$psill[1]
nu.mat <- mat.variog.water$kappa[2]


### Here we compute the weighted sum of squares for the fitted exponential semivariogram
# This is the vector with the weights evaluated at the distances h_m used to compute the empirical semivariogram
weights.exp <- emp.variog.water$np/((exp.variog.water.2)^2)
# This is the vector with the squared difference between the values of the empirical semivariogram at the distances h_m and the
# the values of the fitted exponential semivariogram at the distances h_m
squared.diff.exp <- (emp.variog.water$gamma-exp.variog.water.2)^2
wss.exp <- sum(weights.exp*squared.diff.exp)
wss.exp

