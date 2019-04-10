

######################################################################################
##   This file shows how to fit a Bayesian hierarchical model to univariate
##   spatial data: Y(s) using the spBayes package in R
######################################################################################

# To illustrate how the package works, we will use the dataset BEF.dat
# provided with the package spBayes
# The dataset contains forest inventory data for years 1991 and 2002
# for the Bartlett Experimental Forest located in Bartlett, NH.
# The dataset includes 208 variables among which:
# - specific basal area
# - biomass
# - coordinate of the inventory plot
# - slope
# - elevation
# - tasseled cap brightness (TC1)
# - greenness (TC2)
# - wetness (TC3)
# TC1, TC2, and TC3 are derived from Landsat images collected in Spring, Summer and Fall 2002


## Here I am loading all the packages 
library(spBayes)
library(fields)
library(akima)
library(MBA)
library(MASS)
library(RColorBrewer)
library(classInt)
library(gstat)
library(geoR)


### Reading in the data

data(BEF.dat)

dim(BEF.dat)
n <- dim(BEF.dat)[1]

# To reduce the computational burden for the moment, we just consider a subset of 400 observations from
# the entire dataset.
# We choose 400 observations at random and we consider only the observations for which the variable BE_02BAREA is greater than 0.
# This the selected subset that we are going to consider
set.seed(3028)
index.train <- sample(1:n,400,replace=F)
index.test <- (1:n)[-index.train]
BEF.train <- BEF.dat[index.train,]
index.subset <- which(BEF.train$BE_02BAREA > 0)

BEF.subset <- BEF.train[index.subset,]

# We analyze the variable BE_02BAREA in the BEF.subset on the log scale
Y <- BEF.subset$BE_02BAREA
elev <- BEF.subset$ELEV
spr02.tc1 <- BEF.subset$SPR_02_TC1
spr02.tc2 <- BEF.subset$SPR_02_TC2
spr02.tc3 <- BEF.subset$SPR_02_TC3
x.coord <- BEF.subset$XUTM
y.coord <- BEF.subset$YUTM


# Here we do an exploratory analysis of the data: histogram, Q-Q plot of the data
# and we examine the relationship between the dependent variable and the independent 
# variable elevation, and the three measures of forest canopy given by the satellite 
# (spr02.tc1, spr02.tc2, spr02.tc3)

# The line below is to make multiple plots in one page (the mfrow command)
# and to define the margins of the figure (the mai command)
par(mfrow=c(2,1),mai=c(1,1,0.5,0.5))
hist(Y,breaks=30,col="grey",xlab="Log basal area",prob=TRUE,main="Histogram of 
log basal area",cex.lab=2,cex.axis=2)
qqnorm(Y,main="Q-Q plot of log basal area",cex.lab=2,cex.axis=2,col="red",lwd=2,lty=1)
qqline(Y,col="black")


#
par(mfrow=c(2,2))
plot(elev,Y,type="p",pch=20,col="black",xlab="Elevation",ylab="Log-density",cex.lab=2,cex.axis=2)
plot(spr02.tc1,Y,type="p",pch=20,col="black",xlab="Total canopy 1",ylab="Log-density",cex.lab=2,cex.axis=2)
plot(spr02.tc2,Y,type="p",pch=20,col="black",xlab="Total canopy 2",ylab="Log-density",cex.lab=2,cex.axis=2)
plot(spr02.tc3,Y,type="p",pch=20,col="black",xlab="Total canopy 3",ylab="Log-density",cex.lab=2,cex.axis=2)



# Spatial plot of the observed values and interpolated surfaces
plotvar <- Y
nclr <- 9
plotclr <- brewer.pal(nclr,"YlOrRd")

class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))
colcode <- findColours(class,plotclr)

length.x <- 25
length.y <- 25

x.grid <- seq(min(x.coord),max(x.coord),length=length.x)
y.grid <- seq(min(y.coord),max(y.coord),length=length.y)
z.grid <- matrix(NA,nrow=length.x,ncol=length.y)
z.lim <- range(plotvar)

# Plot of the observed values
par(mfrow=c(1,1),mai=rep(1,4))
image.plot(x.grid,y.grid,z.grid,xlab="X",ylab="Y",zlim=z.lim,
col=plotclr,breaks=class$brks,main="Log density",cex.main=2,cex.axis=2,cex.lab=2)
points(x.coord,y.coord,col=colcode,pch=20,cex=2)


surf.logdens <- interp(x.coord,y.coord,Y)

# Plot of the interpolated surface of log density
image.plot(surf.logdens$x,surf.logdens$y,surf.logdens$z,col=plotclr,
breaks=seq(min(plotvar),max(plotvar),length=nclr+1),
zlim=z.lim,xlab="X",ylab="Y",
main="Interpolated log density",cex.main=2,cex.axis=2,cex.lab=2)
contour(surf.logdens$x,surf.logdens$y,surf.logdens$z,nlevels=10,add=T,col="black")



#####  Empirical spatial EDA

# Fitting a linear model on elevation and the three measures of forest canopy
mod <- lm(Y~elev+spr02.tc1+spr02.tc2+spr02.tc3)
summary(mod)
Call:
lm(formula = Y ~ elev + spr02.tc1 + spr02.tc2 + spr02.tc3)

Residuals:
Min       1Q   Median       3Q      Max 
-2.35080 -0.35326  0.06227  0.50832  1.46133 

Coefficients:
Estimate Std. Error t value Pr(>|t|)    
(Intercept) -0.8410852  2.1148893  -0.398 0.691347    
elev         0.0022888  0.0005762   3.972 0.000105 ***
spr02.tc1   -0.0388430  0.0249854  -1.555 0.121873    
spr02.tc2    0.0695938  0.0158778   4.383 2.03e-05 ***
spr02.tc3   -0.0738726  0.0216058  -3.419 0.000785 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

Residual standard error: 0.6925 on 172 degrees of freedom
Multiple R-squared: 0.4748,	Adjusted R-squared: 0.4626 
F-statistic: 38.88 on 4 and 172 DF,  p-value: < 2.2e-16 


## Plot of the interpolated surface of the residuals to the linear model
res.mod <- mod$residuals
surf.res <- interp(x.coord,y.coord,res.mod)

image.plot(surf.res$x,surf.res$y,surf.res$z,col=heat.colors(100),
zlim=range(surf.res$z,na.rm=T),xlab="X",ylab="Y",
main="Interpolated residuals of log density",cex.main=2,cex.axis=2,cex.lab=2)
contour(surf.res$x,surf.res$y,surf.res$z,nlevels=10,add=T,col="black")


# Computing the empirical semi-variogram for the residuals
logdens.data <- data.frame(cbind(x.coord,y.coord,res.mod))
coordinates(logdens.data)=~x.coord+y.coord

# This is using the default definition of the bins 
emp.variog.logdens.default <- variogram(res.mod~1,logdens.data)

# Plot of the empirical semi-variogram
par(mai=rep(1,4))
plot(emp.variog.logdens.default$dist,emp.variog.logdens.default$gamma,
xlab="Distance",ylab="Semi-variance",main="Empirical semi-variogram of residuals",
cex.lab=2,cex.axis=2,cex.main=2,type="p",col="blue",pch=1,cex=2,ylim=c(0,0.85))


# This is defining the bins by hand: we can do this using the option boundaries in the variogram function
# specifying the bin boundaries explicitly
emp.variog.logdens <- variogram(res.mod~1,logdens.data,boundaries=c(60,seq(150,300,by=50),seq(400,1000,by=100)))


# Plot of the empirical semi-variogram
par(mai=rep(1,4))
plot(emp.variog.logdens$dist,emp.variog.logdens$gamma,
xlab="Distance",ylab="Semi-variance",main="Empirical semi-variogram of residuals",
cex.lab=2,cex.axis=2,cex.main=2,type="p",col="blue",pch=1,cex=2,ylim=c(0,0.85))



# Fitting an exponential semi-variogram to the empirical via Weighted Least Squares
exp.variog.logdens <- fit.variogram(emp.variog.logdens,vgm(psill=0.5,"Exp",range=400,nugget=0.05),fit.method=2)
exp.variog.logdens
model     psill    range
1   Nug 0.0000000  0.00000
2   Exp 0.4349752 35.27601


sigma2.exp <- exp.variog.logdens$psill[2]
tau2.exp <- exp.variog.logdens$psill[1]
phi.exp <- exp.variog.logdens$range[2]

vgm.exp <- variogramLine(vgm(sigma2.exp, "Exp", phi.exp, tau2.exp), dist_vector = seq(0,max(emp.variog.logdens$dist),length=100))

# Plot of the empirical semivariogram and the fitted exponential semi-variogram where the fitting has been done via WLS
par(mai=c(1,1,1,1))
plot(emp.variog.logdens$dist,emp.variog.logdens$gamma,type="p",col="blue",pch=1,
cex=1.5,xlab="Distance",ylab="Semi-variance",
cex.lab=2,cex.main=2,cex.axis=2,ylim=c(0,0.6),
xlim=c(0,max(emp.variog.logdens$dist)))
lines(vgm.exp$dist,vgm.exp$gamma,type="l",col="red",lwd=2,lty=1)



## Here we do maximum likelihood estimation of all the model parameters using Restricted Maximum Likelihood (REML)
## We use as initial values for the covariance parameters the WLS estimates of the parameters 
## Since our model for the mean function uses elevation and the three canopy measures as covariates
## the geodata object that we create needs to contain the coordinate, the dependent variable and the 
## covariates 
logdens.df <- data.frame(cbind(x.coord,y.coord,Y,elev,spr02.tc1,spr02.tc2,spr02.tc3),nrow=length(x.coord),ncol=7)
logdens.geo <- as.geodata(logdens.df,coords.col=c(1,2),data.col=3)

## Fitting the spatial model via REML
logdens.reml <- likfit(logdens.geo, trend = ~elev+spr02.tc1+spr02.tc2+spr02.tc3,cov.model="exponential",
ini=c(sigma2.exp, phi.exp), nugget=tau2.exp, fix.nug = FALSE, lik.met="REML")

likfit: estimated model parameters:
beta0        beta1        beta2        beta3        beta4        tausq      sigmasq          phi 
"   -0.9466" "    0.0016" "   -0.0314" "    0.0601" "   -0.0670" "    0.4259" "    6.9483" "75264.8354" 
Practical Range with cor=0.05 for asymptotic range: 225473.3

likfit: maximised log-likelihood = -178.7

sigma2.reml <- logdens.reml$sigmasq
tau2.reml <- logdens.reml$tausq
phi.reml <- logdens.reml$phi


## Obtaining the REML residuals as a difference of the observed values minus the fitted values
## The fitted values are obtained by taking the estimates of the coefficients times the covariates
res.reml <- Y-as.numeric(logdens.reml$beta)[1]-as.numeric(logdens.reml$beta)[2]*elev-as.numeric(logdens.reml$beta)[3]*spr02.tc1-as.numeric(logdens.reml$beta)[4]*spr02.tc2-as.numeric(logdens.reml$beta)[5]*spr02.tc3

## Here we construct the empirical variogram of the REML residuals 
logdens.reml.data <- data.frame(cbind(x.coord,y.coord,res.reml))
coordinates(logdens.reml.data)=~x.coord+y.coord
emp.variog.logdens.reml <- variogram(res.reml~1,logdens.reml.data,boundaries=c(60,seq(150,300,by=50),seq(400,1000,by=100)))

## This is to compute the exponential variogram with parameters set equal to the REML estimates of the covariance parameters
vgm.exp.reml <- variogramLine(vgm(sigma2.reml, "Exp", phi.reml, tau2.reml), dist_vector = seq(0,max(emp.variog.logdens.reml$dist),length=100))

## Plot of the empirical semivariogram of the REML residuals and the fitted exponential semivariogram with REML estimated parameters
par(mai=c(1,1,1,1))
plot(emp.variog.logdens.reml$dist,emp.variog.logdens.reml$gamma,type="p",col="blue",pch=1,
cex=1.5,xlab="Distance",ylab="Semi-variance",
cex.lab=2,cex.main=2,cex.axis=2,ylim=c(0,0.6),
xlim=c(0,max(emp.variog.logdens.reml$dist)))
lines(vgm.exp.reml$dist,vgm.exp.reml$gamma,type="l",col="red",lwd=2,lty=1)



#########  Spatial analysis from a Bayesian perspective (the code is about the same as the code in the Bayesian hierarchical spatial file

######################################################################################
##   We are analyzing the data using a classical 
##   geostatical Bayesian hierarchical model where we 
##   marginalize the process eta^\prime.
##   The likelihood is then:
##        Y | beta, tau^2, theta ~ MVN(X*beta, Sigma(theta) + tau2 I_n)
##
##   where we model the covariance function using an exponential covariance function
##   Since the exponential covariance function is stationary and isotropic, we can
##   write Sigma(theta) as theta_1 R(theta_2) with theta_1 sill and R(theta_2) 
##   correlation matrix that depends on the distance between points and theta_2
######################################################################################


#############  Priors distribution and hyperparameters  #####################

# theta_1: We specify an Inverse Gamma prior for the sill, theta_1
# with prior mean equal to 30 and infinite variance

# tau^2: We use the same prior for the nugget effect, tau^2


# To determine the prior for theta_2, we look at the distances among points
# This is the matrix with the coordinate of the points
coords <- as.matrix(cbind(BEF.subset$XUTM, BEF.subset$YUTM),nrow=length(BEF.subset$XUTM),ncol=2)

# Since the largest distance among points in the entire dataset is approximately 3656 meters,
# we specify as prior for the theta2 parameter a uniform distribution such that the effective range
# can range between 50 meters to 3700 meters.
# NOTE: In the spBayes, the exponential covariance function is parameterized as: 
#  C(h; theta) = sigma^2*exp(-phi * h)
# and the prior is given in terms of the effective range

# Finally, as prior for beta, we use an improper flat prior


#############  Initial values and tuning parameters for Metropolis-Hastings  #####################

beta.ini <- rep(0,5)
sigma2.ini <- 1/rgamma(1,2,1)
tau2.ini <- 1/rgamma(1,2,0.5)
phi.ini <- runif(1,min=0.0008,max=0.06)


#############  Fitting the model

# This fits the model
model.1 <- spLM(Y~BEF.subset$ELEV+BEF.subset$SPR_02_TC1 + BEF.subset$SPR_02_TC2+ BEF.subset$SPR_02_TC3, 
coords=coords,starting=list("phi"=phi.ini,"sigma.sq"=sigma2.ini, "tau.sq"=tau2.ini,"beta"=beta.ini),
tuning=list("phi"=0.001, "sigma.sq"=0.05, "tau.sq"=0.01),
priors=list("phi.Unif"=c(0.0008, 0.06), "sigma.sq.IG"=c(2, 1),
"tau.sq.IG"=c(2, 0.5),"beta.Flat"), cov.model="exponential",
n.samples=5000, verbose=TRUE, n.report=100)

# To see the output generated by the function spLM, we can use
names(model.1)
model.1$p.theta.samples[1:5,]

# This produces traces plots and plots of the posterior marginal distribution of the covariance parameters
par(mai=rep(0.4,4))
plot(model.1$p.theta.samples[,1:3])

n.samples <- 5000
burn.in <- 0.5*n.samples
model.1.other.pars <- spRecover(model.1, start=burn.in, verbose=FALSE)


names(model.1.other.pars)
[1] "p.theta.samples"         "acceptance"              "tuning"                  "Y"                      
[5] "X"                       "coords"                  "is.pp"                   "cov.model"              
[9] "nugget"                  "beta.prior"              "beta.Norm"               "x.names"                
[13] "run.time"                "p.beta.recover.samples"  "p.theta.recover.samples" "p.w.recover.samples"    

# This produces traces plots and plots of the posterior marginal distribution of the beta coefficients
dim(model.1.other.pars$p.beta.recover.samples)

par(mai=rep(0.4,4),mfrow=c(2,2))
plot(model.1.other.pars$p.beta.recover.samples[,1:4])

par(mai=rep(0.4,4))
plot(model.1.other.pars$p.beta.recover.samples[,5])

par(mai=rep(0.4,4),mfrow=c(2,2))
plot(model.1.other.pars$p.theta.samples[,1:3])


# To obtain estimates and 95% credible of the covariance parameters
round(summary(model.1.other.pars$p.theta.samples)$quantiles[,c(1,3,5)],2)

# To obtain estimates and 95% credible of the beta coefficients
round(summary(model.1.other.pars$p.beta.recover.samples)$quantiles[,c(1,3,5)],2)



# To obtain the spatial process
eta.summary <- summary(mcmc(t(model.1.other.pars$p.w.recover.samples)))$quantiles[,c(1,3,5)]
eta.summary


surf.eta <- interp(x.coord,y.coord,as.numeric(eta.summary[,1]))

par(mfrow=c(1,1))
image.plot(surf.eta$x,surf.eta$y,surf.eta$z,col=heat.colors(100),
zlim=range(surf.eta$z,na.rm=T),xlab="X",ylab="Y",
main="Interpolated posterior mean of spatial process",cex.main=2,cex.axis=2,cex.lab=2)
contour(surf.eta$x,surf.eta$y,surf.eta$z,nlevels=10,add=T,col="black")


#############  Making prediction

# Here we want to make prediction at sites without observations
# We make predictions at the remaining 237 sites

n.pred <- 237

# This is the subset of the BEF.dat dataset that does not contain the data already used to fit the model
BEF.pred <- BEF.dat[c(201:n),]

# In this matrix, we place the covariates for the sites where we want to make predictions
BEF.predcov <- matrix(cbind(rep(1, n.pred), BEF.pred$ELEV,BEF.pred$SPR_02_TC1, BEF.pred$SPR_02_TC2,BEF.pred$SPR_02_TC3),nrow=n.pred,ncol=5)

# These are the coordinates of the sites where we want to make prediction
pred.coords <- cbind(BEF.pred$XUTM, BEF.pred$YUTM)

# This command generates predictions at the 237 selected sites for prediction, by using composition
# sampling and the MCMC samples generated earlier (here we are prediciting from the posterior predictive distribution). 
# Since we use 1000 as burnin, we are starting to make predictions from the 2500-th iterations
# and we make predictions every 2 iterations (this is the thinning parameter)
pred <- spPredict(model.1, pred.coords=pred.coords, pred.covars=BEF.predcov, start=burn.in, thin=2)
names(pred)

dim(pred$p.y.predictive.samples)

## Here we compute the posterior mean of the predictions. 
post.pred.mean <- rowMeans(pred$p.y.predictive.samples)
post.pred.mean[1:10]

## Here we compute the 95% posterior predictive intervals at the 237 sites
post.pred.95ci <- apply(pred$p.y.predictive.samples,1,quantile,c(0.025,0.975))
post.pred.95ci[,1:10]


Y.pred <- log(BEF.pred$BE_02BAREA)

### Assessing the quality of the predictions: Mean absolute error
mean(abs((post.pred.mean-Y.pred))[Y.pred!=-Inf])
[1] 0.746323

## Root mean squared error
sqrt(mean(((post.pred.mean-Y.pred)^2)[Y.pred!=-Inf]))
[1] 0.9502872

