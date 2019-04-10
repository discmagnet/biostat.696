

####################################################################################
#####
#####  This code illustrates how to estimate parameters for spatial processes in R
#####  using Maximum Likelihood and Restricted Maximum Likelihood
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

# These are the fitted values and the residuals
mod.water <- lm(water.soil~x.soil+y.soil+x.soil.sq+y.soil.sq+x.y.soil+x.soil.cub+x.soil.sq.y.soil+x.soil.y.soil.sq+y.soil.cub)
trend.water <- mod.water$fitted.values
res.water <- water.soil-trend.water

# Here we are computing the empirical variogram using the variogram function in gstat
water.data <- data.frame(cbind(x.soil,y.soil,res.water))
emp.variog.water <- variogram(res.water~1,locations=~x.soil+y.soil,data=water.data)


# Here we fit the spherical semi-variogram to the residuals of the pH of water
sph.variog.water <- fit.variogram(emp.variog.water,vgm(psill=0.01,"Sph",range=10,nugget=0.01),fit.method=2)
sph.variog.water

sigma2.sph <- sph.variog.water$psill[2]
phi.sph <- sph.variog.water$range[2]
tau2.sph <- sph.variog.water$psill[1]
dist.vec <- emp.variog.water$dist


############################################################
#
#     ML and REML estimation of the covariance function paraeters
#
############################################################

## To fit a spatial model using likelihood based methods
## we use the likfit function in the geoR library

library(geoR)

## In order to use the likfit function we need to put the data in the geodata format
## Thus, we first create a data frame and then we put it in the geodata format
## by specifying which ones are the coordinate columns
water.df <- data.frame(cbind(x.soil,y.soil,water.soil),nrow=length(x.soil),ncol=3)
water.geo <- as.geodata(water.df,coords.col=c(1,2),data.col=3)

## To use the likfit function, we need to specify initial values for the parameters
## As initial values, we use the estimates of the spherical semi-variogram parameters obtained via WLS
sigma2.sph <- sph.variog.water$psill[2]
phi.sph <- sph.variog.water$range[2]
tau2.sph <- sph.variog.water$psill[1]
dist.vec <- emp.variog.water$dist

## In the likfit (type ?likfit for more details), we can specify the model for the spatial trend,
## we need to provide initial estimates for the covariance parameters
## and we can specify whether we want to perform MLE or REML estimation (in the lik.met entry) 

## This is using ML
water.ml <- likfit(water.geo, trend = ~x.soil+y.soil+I(x.soil^2)+I(y.soil^2)+I(x.soil*y.soil)+I(x.soil^3)+
I(x.soil^2*y.soil)+I(x.soil*y.soil^2)+I(y.soil^3), 
cov.model="spherical",ini=c(sigma2.sph, phi.sph), nugget=tau2.sph, fix.nug = FALSE, lik.met="ML")

# Here we look at the output
water.ml
## These are the estimates of the trend parameters (beta) and of the covariance parameters (nugget effect, tau^2, partial sill, theta1, 
## and theta2) using MLE
#likfit: estimated model parameters:
#beta0     beta1     beta2     beta3     beta4     beta5     beta6     beta7     beta8     beta9     tausq   sigmasq
#" 5.7608" " 0.0194" " 0.0102" "-0.0014" "-0.0001" "-0.0002" " 0.0000" " 0.0000" " 0.0000" " 0.0000" " 0.0035" " 0.0200" 
#phi 
#"19.5330" 
#Practical Range with cor=0.05 for asymptotic range: 19.53302

#likfit: maximised log-likelihood = 182.8


## This is using REML
water.reml <- likfit(water.geo, trend = ~x.soil+y.soil+I(x.soil^2)+I(y.soil^2)+I(x.soil*y.soil)+I(x.soil^3)+
I(x.soil^2*y.soil)+I(x.soil*y.soil^2)+I(y.soil^3), 
cov.model="spherical",ini=c(sigma2.sph, phi.sph), nugget=tau2.sph, fix.nug = FALSE, lik.met="REML")

# Here we look at the output
water.reml

## These are the estimates of the trend parameters (beta) and of the covariance parameters (nugget effect, tau^2, partial sill, theta1, 
## and theta2) using REML
#likfit: estimated model parameters:
#beta0     beta1     beta2     beta3     beta4     beta5     beta6     beta7     beta8     beta9     tausq   sigmasq 
#" 5.7628" " 0.0178" " 0.0102" "-0.0013" "-0.0001" "-0.0003" " 0.0000" " 0.0000" " 0.0000" " 0.0000" " 0.0014" " 0.0292" 
#phi 
#"20.6379" 
#Practical Range with cor=0.05 for asymptotic range: 20.63786

#likfit: maximised log-likelihood = 181.3


############################################################
#
#     AIC/BIC
#
############################################################


## This is to test whether we should model the covariance matrix using a spherical covariance 
## function or using a matern covariance function.
## In both models, the spatial trend is a cubic polynomial.
## Since these models are not nested, we cannot use the likelihood ratio test to determine
## which model is more appropriate to the data, but we can compare the AIC and BIC for the two models
## The model with the lower AIC (and BIC) is the best-fitting model.

water.reml.sph <- likfit(water.geo, trend = ~x.soil+y.soil+I(x.soil^2)+I(y.soil^2)+I(x.soil*y.soil)+I(x.soil^3)+
I(x.soil^2*y.soil)+I(x.soil*y.soil^2)+I(y.soil^3), 
cov.model="spherical",ini=c(sigma2.sph, phi.sph), nugget=tau2.sph, fix.nug = FALSE, lik.met="REML")

names(water.reml.sph)

# The likfit returns an object that also has the AIC and BIC value for the model. 
# To access them, simply subset in the output list created by the likfit function
# AIC and BIC values corresponding to the fitted model
water.reml.sph$AIC
water.reml.sph$BIC

#


# Here we fit the Matern semi-variogram to the residuals of the pH of water
mat.variog.water <- fit.variogram(emp.variog.water,vgm(psill=0.01,"Mat",range=10,nugget=0.01,kappa=1),fit.method=2)
mat.variog.water

sigma2.mat <- mat.variog.water$psill[2]
phi.mat <- mat.variog.water$range[2]
tau2.mat <- mat.variog.water$psill[1]
nu.mat <- mat.variog.water$kappa[2]


water.reml.mat <- likfit(water.geo, trend = ~x.soil+y.soil+I(x.soil^2)+I(y.soil^2)+I(x.soil*y.soil)+I(x.soil^3)+
I(x.soil^2*y.soil)+I(x.soil*y.soil^2)+I(y.soil^3),  
cov.model="matern",ini=c(sigma2.mat, phi.mat), nugget=tau2.mat, kappa=nu.mat,fix.nug = FALSE, lik.met="REML")

water.reml.mat$AIC
water.reml.mat$BIC


