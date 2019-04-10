

####################################################################################
#####
#####  This code illustrates how to simulate Gaussian spatial processes with
#####  a given covariance function in R
#####
####################################################################################


## Here we first make plots of some parametric covariance functions.

## The function "show.vgm" in the gstat library makes plots of different parametric 
## semi-variogram model. For more information, first load the gstat library and then type 
## ?show.vgm

## To simulate Gaussian spatial processes with a given covariance function in R, we need to install and upload the MASS library

## Here we set up the coordinates of the points where we want to simulate the spatial process
## We use a 21x21 grid of points  

## These are the unique set of x and y coordinates 
x.grid <- seq(0,1,by=0.05)
y.grid <- seq(0,1,by=0.05)
n.grid <- length(x.grid)

## These are the coordinates of the points: we set up a matrix with 21*21
## rows and 2 columns. Each row gives the coordinate of each point in the grid
coord.pts <- matrix(0,nrow=n.grid*n.grid,ncol=2)

## Here we fill the matrix of coordinates by column
coord.pts[,1] <- rep(x.grid,each=n.grid)
coord.pts[,2] <- rep(y.grid,n.grid)


## This is the function that computes the covariance
## between two points when the covariance function of the spatial
## process is modeled as a stationary isotropic:
## 1. Gaussian covariance function
## 2. Exponential covariance function

gauss.cov <- function(sigma2,phi,tau2,dist.pts){
	
	cov.value <- tau2+(sigma2*exp(-(dist.pts^2)/(phi^2)))
	return(cov.value)
		
}

exp.cov <- function(sigma2,phi,tau2,dist.pts){
	
	cov.value <- tau2+(sigma2*(exp(-(dist.pts/phi))))
	return(cov.value)
	
}

## Note that there are also built-in functions in the fields package or in the gstat that can compute these covariance functions.
## See for example: Exp.cov or stationary.cov for more details in the fields package
## or vgm in the gstat package

## Here we create the matrix with distances among points.
## This is a N by N matrix where the (i,j) component is the distance between the i-th point s_i
## and the j-th point s_j

N <- dim(coord.pts)[1]
## This is doing it by hand
Dist.mat <- matrix(0,nrow=N,ncol=N)

for(i in 1:N){
	for(j in 1:N){
		Dist.mat[i,j] <- sqrt(sum((coord.pts[i,]-coord.pts[j,])^2))	
	}
}


## This is using the rdist.function in the fields package
library(fields)
Dist.mat <- rdist(coord.pts,coord.pts)


## Here we construct the N by N covariance matrix whose (i,j) element is the 
## covariance between Y(s_i) and Y(s_j).
## 1. Sigma.gauss is the N by N covariance matrix when we model the covariance 
## function for the spatial process using a spherical covariance function
##
## 2. Sigma.exp is the N by N covariance matrix when we model the covariance 
## function for the spatial process using the exponential covariance function

Sigma.gauss <- matrix(0,nrow=N,ncol=N)
Sigma.exp <- matrix(0,nrow=N,ncol=N)

## Here we specify the values for the parameters of the two covariance function
sigma2.gauss <- 1
phi.gauss <- 3
tau2.gauss <- 0

sigma2.exp <- 1
phi.exp <- 2
tau2.exp <- 0

for(i in 1:N){
	for(j in 1:N){
		Sigma.gauss[i,j] <- gauss.cov(sigma2.gauss,phi.gauss,tau2.gauss,Dist.mat[i,j])
		Sigma.exp[i,j] <- exp.cov(sigma2.exp,phi.exp,tau2.exp,Dist.mat[i,j])
	}
}

Sigma.gauss[1:3,1:3]
Sigma.exp[1:3,1:3]



## Here we simulate realizations of two Gaussian spatial processes with mean 0
## and covariance functions given by the Gaussian covariance function
## and the exponential covariance function with the parameters specified above
## To do so, we use the mvrnorm function in the MASS package

library(MASS)

field.gauss <- mvrnorm(1,rep(0,N),Sigma.gauss)
field.exp <- mvrnorm(1,rep(0,N),Sigma.exp)

## Here we plot the two simulated Gaussian spatial processes
image.plot(x.grid,y.grid,matrix(field.gauss,n.grid,n.grid),col=terrain.colors(100),xlab="X",ylab="Y",main="Gaussian covariance function")
image.plot(x.grid,y.grid,matrix(field.exp,n.grid,n.grid),col=terrain.colors(100),xlab="X",ylab="Y",main="Exponential covariance function")


########

N <- dim(coord.pts)[1]

## Here we simulate a realization of a Gaussian spatial processes with mean 0
## and covariance functions given by a geometrically anisotropic exponential covariance function 
## corresponding to a rotation of 30 degree 

A <- matrix(c(sqrt(3)/2,-1/2,1/2,sqrt(3)/2),nrow=2,ncol=2,byrow=T)

## Here we derive the matrix h^prime*A*h: these represents the distances in the new set of coordinates (the coordinate system obtained 
## by rotating the old ones by 30 degrees)
Dist.geo.aniso <- matrix(0,nrow=N,ncol=N)
for(i in 1:N){
	for(j in 1:N){
		# this is the vector h corresponding to the separation between point s_i and point s_j 
		h.vec <- c(coord.pts[i,1]-coord.pts[j,1],coord.pts[i,2]-coord.pts[j,2])
		# here we compute (h^prime*A*h)^1/2
		Dist.geo.aniso[i,j] <- sqrt(t(matrix(h.vec,nrow=2,ncol=1))%*%A%*%matrix(h.vec,nrow=2,ncol=1))
	}
}


Sigma.geo.aniso.exp <- matrix(0,nrow=N,ncol=N)

for(i in 1:N){
	for(j in 1:N){
		Sigma.geo.aniso.exp[i,j] <- exp.cov(sigma2.exp,phi.exp,tau2.exp,Dist.geo.aniso[i,j])
	}
}

Sigma.geo.aniso.exp[1:3,1:3]

field.geo.aniso.exp <- mvrnorm(1,rep(0,N),Sigma.geo.aniso.exp)

## Here we plot the simulated Gaussian spatial processes with geometrically anisotropic exponential covariance function
image.plot(x.grid,y.grid,matrix(field.geo.aniso.exp,n.grid,n.grid),col=terrain.colors(100),xlab="X",ylab="Y",main="Geometrically Anisotropic \n Exponential covariance function")
