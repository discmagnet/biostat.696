

####################################################################################
#####
#####  This code illustrates how to perform kriging in R
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


#########################################################################
#########################   Simple kriging  #############################
#########################################################################
# We want to predict the pH of water at the site s0=(x0,y0)=(13,13)
# We assume that the mean mu(s) of the spatial process Y(s) is known and equal to 5.8.
# We also assume that the small-scale variation is second-order stationary with 
# spherical covariance function with parameters sigma^2=0.02, phi=20.0 and tau2=0.001

x0 <- 13
y0 <- 13

###
# Here we derive the simple kriging predictor by setting up the equations 

# This is the c vector that contains the covariance between Y(s0) and Y(si) with i=1,..,250
c.vec <- rep(0,length(x.soil))

# Here we compute the distances between the sites where we have observations
Dist.mat <- rdist(matrix(cbind(x.soil,y.soil),nrow=length(x.soil),ncol=2))

# This instead is the vector that contains the distance between s0 and points si, i=1,..,250
dist.pt <- rdist(matrix(c(x0,y0),nrow=1,ncol=2),matrix(cbind(x.soil,y.soil),nrow=length(x.soil),ncol=2))

# This is the function that computes the covariance function
spher.cov <- function(sigma2,phi,tau2,dist){
	cov.out <- NULL
	if(dist >= phi){
		cov.out <- 0
	}
	if(dist < phi){
		cov.out <- sigma2*(1-(((3*dist)/(2*phi))-((dist^3)/(2*(phi^3)))))
	}
	if(dist==0){
		cov.out <- tau2+sigma2	
	}
	return(cov.out)
}


# Here we fill the c vector 
for(i in 1:length(c.vec)){
	c.vec[i] <- spher.cov(0.02,20.0,0.001,dist.pt[i])
}

# Here we fill the Sigma matrix
Sigma <- matrix(0,length(x.soil),length(y.soil))
for(i in 1:length(x.soil)){
	for(j in 1:length(y.soil)){
		Sigma[i,j] <- spher.cov(0.02,20.0,0.001,Dist.mat[i,j])
	}
}


# We assume that the mean is known and is contant and equal to 5.8 everywhere.
# Here we are setting up the mu and the Y vector and we are specifying the value of mu at s0
mu <- rep(5.8,length(water.soil))
mu.s0 <- 5.8
Y <- water.soil

# This is deriving the simple kriging predictor at s0 using the simple kriging formula: p(Y; s0)= mu(s0) + t(c)%*%solve(Sigma)%*%(Y-mu)
# with Y vector with all the 250 observations and mu vector with the mean function at sites si, i=1,..,250
pred.s0 <- mu.s0+t(matrix(c.vec,nrow=length(x.soil),ncol=1))%*%(solve(Sigma))%*%(matrix(Y-mu,nrow=length(Y),ncol=1))
pred.s0
[1] 5.865278

sk.variance.s0 <- (0.02+0.001)-(t(matrix(c.vec,nrow=length(x.soil),ncol=1))%*%(solve(Sigma))%*%matrix(c.vec,nrow=length(x.soil),ncol=1))
sk.variance.s0
[1,] 0.005023182


### 
# Here we derive the simple kriging predictor using geoR

library(geoR)

## As usual, we first create a data frame and then we put it in the geodata format
## by specifying which ones are the coordinate columns
water.df <- data.frame(cbind(x.soil,y.soil,water.soil),nrow=length(x.soil),ncol=3)
water.geo <- as.geodata(water.df,coords.col=c(1,2),data.col=3)

# This specifies the options for the type of kriging that we want to perform
kc.sk.control <- krige.control(type.krige="sk",trend.d="cte",trend.l="cte",beta=mu.s0,cov.model="spherical",
cov.pars=c(0.02,20.0),nugget=0.001)

# Here we putting the coordinates of the point for which we want to obtain the simple kriging predictor in a matrix form:
loc.sk <- matrix(c(x0,y0),nrow=1,ncol=2)
kc.sk.s0 <- krige.conv(water.geo,locations=loc.sk,krige=kc.sk.control)
kc.sk.s0$predict
data 
5.865652 

kc.sk.s0$krige.var
[1] 0.005330207





#########################################################################
#########################   Ordinary kriging  #############################
#########################################################################
# We want to predict the pH of water at the site s0=(x0,y0)=(13,13)
# We assume that the mean mu(s) of the spatial process Y(s) is constant but unknown.
# We also assume that the small-scale variation is second-order stationary with 
# spherical covariance function with parameters that we estimate via WLS


water.data <- data.frame(cbind(x.soil,y.soil,water.soil))
coordinates(water.data)=~x.soil+y.soil

# Note that here since we are assuming that the mean of Y(s) is constant we can just compute 
# the empirical variogram on the data directly
emp.variog.water <- variogram(water.soil~1,water.data)

# Here we fit a spherical variogram to the empirical variogram via WLS
sph.variog.water <- fit.variogram(emp.variog.water,vgm(psill=0.01,"Sph",range=10,nugget=0.01),fit.method=2)
sph.variog.water

# Here we take the estimate of the WLS procedure and we will use them to set up the ordinary kriging equations
sigma2.wls <- sph.variog.water$psill[2]
phi.wls <- sph.variog.water$range[2]
tau2.wls <- sph.variog.water$psill[1]


###
# Here we derive the ordinary kriging predictor by setting up the equations 
# This is the vector of all 1s (it is specified as a matrix for computation purposes in R)
one.matrix <- matrix(rep(1,length(x.soil)),nrow=length(x.soil),ncol=1)

# This is the vector that contains the covariance between s0 and the other points
c.vec.wls <- rep(0,length(c.vec))
for(i in 1:length(c.vec)){
	c.vec.wls[i] <- spher.cov(sigma2.wls,phi.wls,tau2.wls,dist.pt[i])
}

# This is the covariance matrix Sigma between pH of water at the 250 sites si
Sigma.wls <- matrix(0,length(x.soil),length(y.soil))
for(i in 1:length(x.soil)){
	for(j in 1:length(y.soil)){
		Sigma.wls[i,j] <- spher.cov(sigma2.wls,phi.wls,tau2.wls,Dist.mat[i,j])
	}
}

# Here we obtain the lambda vector using the formula: 
#      lambda= t(c + one.matrix%*%(1-t(one.matrix)%*%solve(Sigma)%*%c)%*%solve(t(one.matrix)%*%solve(Sigma)%*%one.matrix))%*%solve(Sigma)
lambda.wls <- t(matrix(c.vec.wls,nrow=length(x.soil),ncol=1)+(one.matrix%*%(1-t(one.matrix)%*%solve(Sigma.wls)%*%matrix(c.vec.wls,nrow=length(x.soil),ncol=1))%*%solve(t(one.matrix)%*%solve(Sigma.wls)%*%one.matrix)))%*%solve(Sigma.wls)

# Here we compute the ordinary kriging predictor at s0 as the sum of the observations each weighted using the lambda weights
pred.s0.ok.wls <- sum(lambda.wls*Y)
pred.s0.ok.wls
[1] 5.877536

# Here we compute the ordinary kriging prediction variance
part1.ok.var.wls <- (sigma2.wls+tau2.wls)-(matrix(c.vec.wls,nrow=1,ncol=length(x.soil))%*%(solve(Sigma.wls))%*%t(matrix(c.vec.wls,nrow=1,ncol=length(x.soil))))
part2.ok.var.wls <- t((1-t(one.matrix)%*%solve(Sigma.wls)%*%matrix(c.vec.wls,nrow=length(x.soil),ncol=1)))%*%solve(t(one.matrix)%*%solve(Sigma.wls)%*%one.matrix)%*%(1-t(one.matrix)%*%solve(Sigma.wls)%*%matrix(c.vec.wls,nrow=length(x.soil),ncol=1))

ok.variance.s0 <- part1.ok.var.wls+part2.ok.var.wls
ok.variance.s0
[1,] 0.006427987


### 
# Here we derive the ordinary kriging predictor using geoR

## As usual, we first create a data frame and then we put it in the geodata format
## by specifying which ones are the coordinate columns
water.df <- data.frame(cbind(x.soil,y.soil,water.soil),nrow=length(x.soil),ncol=3)
water.geo <- as.geodata(water.df,coords.col=c(1,2),data.col=3)

# This specifies the options for the type of kriging that we want to perform
kc.ok.control <- krige.control(type.krige="ok",trend.d="cte",trend.l="cte",cov.model="spherical",
cov.pars=c(sigma2.wls,phi.wls),nugget=tau2.wls)

loc.ok <- matrix(c(x0,y0),nrow=1,ncol=2)
kc.ok.s0 <- krige.conv(water.geo,locations=loc.ok,krige=kc.ok.control)
kc.ok.s0$predict
data 
5.877536 

kc.ok.s0$krige.var
[1] 0.006427987


#########################################################################
#########################   Universal kriging  #############################
#########################################################################
# We want to predict the pH of water at the site s0=(x0,y0)=(13,13)
# We assume that the mean mu(s) of the spatial process Y(s) is a linear function 
# of the x and y coordinate, that also contains an intercept term.
# We also assume that the small-scale variation is second-order stationary with 
# spherical covariance function with parameters that we estimate via REML.



### 
# Here we derive the universal kriging predictor using geoR

# First we estimate the covariance parameters using REML assuming that the large-scale variation
# is a linear function of x and y

# As usual, we first fit a linear model to the data, get the residuals and build the empirical variogram
# for the estimated residuals
mod.water.lin <- lm(water.soil~x.soil+y.soil)
trend.water.lin <- mod.water.lin$fitted.values
res.water.lin <- water.soil-trend.water.lin

# Here we derive the WLS estimate of the spherical semi-variogram
water.res.data <- data.frame(cbind(x.soil,y.soil,res.water.lin))
coordinates(water.res.data)=~x.soil+y.soil
emp.variog.res.water <- variogram(res.water.lin~1,water.res.data)
sph.variog.res.water <- fit.variogram(emp.variog.res.water,vgm(psill=0.01,"Sph",range=10,nugget=0.01),fit.method=2)

# There are the WLS estimates of the parameters: we will use these as initial values for the REML procedure
sigma2.sph <- sph.variog.res.water$psill[2]
phi.sph <- sph.variog.res.water$range[2]
tau2.sph <- sph.variog.res.water$psill[1]

# Here we get the REML estimates of the spherical covariance parameters
water.reml.sph <- likfit(water.geo, trend = ~x.soil+y.soil, 
cov.model="spherical",ini=c(sigma2.sph, phi.sph), nugget=tau2.sph, fix.nug = FALSE, lik.met="REML")

water.reml.sph
sigma2.reml <- water.reml.sph$sigmasq
phi.reml <- water.reml.sph$phi
tau2.reml <- water.reml.sph$tausq

# Here we define the control parameters for universal kriging.
# Since we are assuming that the mean is a linear polynomial in the coordinate, we now change the specification of trend.d and trend.l
# and we specify that they have to be a linear polynomial
kc.uk.control <- krige.control(type.krige="ok",trend.d="1st",trend.l="1st",cov.model="spherical",
cov.pars=c(sigma2.reml,phi.reml),nugget=tau2.reml)

loc.uk <- matrix(c(x0,y0),nrow=1,ncol=2)
kc.uk.s0 <- krige.conv(water.geo,locations=loc.uk,krige=kc.uk.control)
kc.uk.s0$predict
data 
5.865465 

kc.uk.s0$krige.var
[1] 0.007055347









