

###############################################################################################
#####
#####
#####           This code shows how to fit a generalized spatial linear model to spatial
#####           data using the spBayes package (which accomodates only Bernouilly/Binomial and
#####           Poisson data as non-normal data).
#####           We will use the spGLM function in the spBayes package.
#####
###############################################################################################


library(RColorBrewer)
library(classInt)
library(geoR)
library(spBayes)
library(fields)
library(MBA)
library(geoRglm)
library(akima)


##################################################################
##################################################################
##
##                    SPATIAL LOGISTIC REGRESSION
##
## Here we are reading in the dataset of malaria prevalence in Gambia (loaded through
## the geoR package).
## The dataset contains information on the prevalence of malaria in villages in Gambia.
## The variables in the dataset are:
##   - x-coordinate of the village
##   - y-coordinate of the village
##   - the binary outcome variable denoting whether a test on a 5-year old child
##     resulted positive (1) or not
##     (0).
##   - the age of the child in days
##   - an indicator (0/1) to denote whether bed nets were used
##   - a satellite-derived measure of the green-ness of the vegetation surronding the village
##   - a variable denoting whether there is (1) or not (0) a health center in the village.
##
##
##################################################################
##################################################################


# Reading the data about malaria in Gambia
data(gambia)
gambia[1:3,]

# This to make a map of the prevalence of malaria in Gambia
# These are the Easting and Northing coordinate of each village.
x.gambia <- gambia$x
y.gambia <- gambia$y

## Changing the UTM coordinates in km (we will use these coordinates in the Bayesian spatial generalized linear mixed model)
x.gambia.km <- x.gambia/1000
y.gambia.km <- y.gambia/1000


# Since the data is at the child level (0 if the child tested negative, 1 otherwise)
# and each child is assigned the village coordinate, we want to get only the village
# geographica coordinates and map the prevalence of malaria in each village

# This puts all the geographical coordinates together (in meters and km)
coord <- cbind(x.gambia,y.gambia)
coord.km <- cbind(x.gambia.km,y.gambia.km)
# This selects only the unique coordinates and thus returns the coordinates of each village
un.coord <- unique(coord)
un.coord.km <- unique(coord.km)

# Now we count how many children were tested in each village and how many cases of malaria
# there were in each village
no.children.village <- rep(0,dim(un.coord)[1])
no.malaria.village <- rep(0,dim(un.coord)[1])

for(i in 1:dim(un.coord)[1]){
    ind.village <- which(coord[,1]==un.coord[i,1] & coord[,2]==un.coord[i,2])
    no.children.village[i] <- length(ind.village)
    no.malaria.village[i] <- length(which(gambia$pos[ind.village]==1))
}

# This derives the prevalence of malaria
prev.gambia <- no.malaria.village/no.children.village

# Here we plot the prevalence of malaria
plotvar <- prev.gambia
nclr <- 9
plotclr <- brewer.pal(nclr,"BuPu")

class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))
colcode <- findColours(class,plotclr)


length.x <- 25
length.y <- 25

x.grid <- seq(min(un.coord[,1]),max(un.coord[,1]),length=length.x)
y.grid <- seq(min(un.coord[,2]),max(un.coord[,2]),length=length.y)
z.grid <- matrix(NA,nrow=length.x,ncol=length.y)
z.lim <- range(plotvar)


par(mfrow=c(1,1),mai=rep(1,4))
image.plot(x.grid,y.grid,z.grid,xlab="X",ylab="Y",zlim=z.lim,
col=plotclr,breaks=class$brks,main="Prevalence of malaria in Gambia",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
points(un.coord[,1],un.coord[,2],col=colcode,pch=20,cex=2)


## Other information in the dataset is the green-ness level in the village
## and whether the child used a bed net or not
## We want to summarize this information to the village level

## We define two new variables, the green-ness at the village and the percentage of
## children that use bed nets in each village
green.coord.gambia <- rep(0,dim(un.coord)[1])
net.use.coord.gambia <- rep(0,dim(un.coord)[1])


for(i in 1:dim(un.coord)[1]){
    
    ## this tells R which observations in the dataset refer to the i-th village
    ind.village <- which(coord[,1]==un.coord[i,1] & coord[,2]==un.coord[i,2])
    ## here we compute the average green-ness level for the i-th village
    green.coord.gambia[i] <- mean(gambia$green[ind.village])
    ## here we compute the percentage of children that use bed nets in the i-th village
    net.use.coord.gambia[i] <- sum(gambia$netuse[ind.village])/length(ind.village)
    print(i)
}

## Here we plot the two variables: green-ness and percentage of children using bed nets
plotvar <- green.coord.gambia
nclr <- 9
plotclr <- brewer.pal(nclr,"Greens")

class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))
colcode <- findColours(class,plotclr)


length.x <- 25
length.y <- 25

x.grid <- seq(min(un.coord[,1]),max(un.coord[,1]),length=length.x)
y.grid <- seq(min(un.coord[,2]),max(un.coord[,2]),length=length.y)
z.grid <- matrix(NA,nrow=length.x,ncol=length.y)
z.lim <- range(plotvar)


par(mfrow=c(1,1),mai=rep(1,4))
image.plot(x.grid,y.grid,z.grid,xlab="X",ylab="Y",zlim=z.lim,
col=plotclr,breaks=class$brks,main="Green-ness level in Gambia",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
points(un.coord[,1],un.coord[,2],col=colcode,pch=20,cex=2)


plotvar <- net.use.coord.gambia
nclr <- 9
plotclr <- brewer.pal(nclr,"Oranges")

class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))
colcode <- findColours(class,plotclr)


length.x <- 25
length.y <- 25

x.grid <- seq(min(un.coord[,1]),max(un.coord[,1]),length=length.x)
y.grid <- seq(min(un.coord[,2]),max(un.coord[,2]),length=length.y)
z.grid <- matrix(NA,nrow=length.x,ncol=length.y)
z.lim <- range(plotvar)


par(mfrow=c(1,1),mai=rep(1,4))
image.plot(x.grid,y.grid,z.grid,xlab="X",ylab="Y",zlim=z.lim,
col=plotclr,breaks=class$brks,main="Prevalence of net use in Gambia",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
points(un.coord[,1],un.coord[,2],col=colcode,pch=20,cex=2)

######################
######################
######################  Fitting a generalized spatial linear model (Binomial/Bernouilli data).
######################  To fit this model, if the data is in the binomial form, we need
######################  to tell R how many trials (called weights) we have
######################  and how many successes at each location.
######################  The trials here are the number of children tested in each village,
######################  while the successes are the number of children who tested positive.
######################  After we have derived those, we will fit a spatial logistic
######################  regression model with covariates, green-ness level in a village
######################  and percentage of bed net use in a village, and with spatial
######################  random effects.
######################  The spatial random effects will be modeled to have an exponential ######################  covariance function without nugget effect.
######################  We will then need to specify priors for the logistic
######################  regression coefficients and for the covariance parameters sigma^2
######################  and phi.
######################  For the logistic regression coefficients (the beta's) we will use
######################  with mean 0 and a very large variance.
######################  For sigma^2, we will use an Inverse Gamma distribution with
######################  shape parameter 2 and scale parameter 1.
######################  For the range, we will use a Uniform distribution such that the
######################  maximum distance between villages is contained in it.
######################

set.seed(3271)

## This is to define the weights in the logistic regression and this is to define
## the number of successes
weights.obs <- no.children.village
y <- no.malaria.village


## As initial values for the logistic regression coefficients, we get the estimated
## regression coefficients in a logistic regression without spatial correlation
logist.reg <- glm((y/weights.obs)~net.use.coord.gambia+green.coord.gambia, weights=weights.obs, family="binomial")
beta.starting <- coefficients(logist.reg)
beta.tuning <- 1.3*t(chol(vcov(logist.reg)))

## This computes the maximum distance between villages in km
max.dist <- max(rdist(un.coord.km))
max.dist

## Remember that in spBayes, instead of specifying a prior for the range, we need
## to specify a prior for 1/range. So, we will use a Uniform distribution such that
## the inverse of 3 times the maximum distance is contained in it (this will mean that the correlation is equal to 0.05 at tha maximum distance).
phi.starting <- runif(1,0.001,1)

## This defines how many iterations we are going to run: 30,000 iterations
## divided in 600 batches of length 50
n.batch <- 1000
batch.length <- 50
n.samples <- n.batch*batch.length


spatial.logist <- spGLM(y~net.use.coord.gambia+green.coord.gambia, family="binomial", coords=un.coord.km,weights=weights.obs,
starting=list("beta"=beta.starting, "phi"=phi.starting,"sigma.sq"=1, "w"=0),
tuning=list("beta"=beta.tuning, "phi"=0.05, "sigma.sq"=0.5, "w"=0.5),
priors=list("beta.Normal"=list(rep(0,3),rep(100,3)), "phi.Unif"=c(0.001,1), "sigma.sq.IG"=c(2, 1)),
amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
cov.model="exponential", verbose=TRUE, n.report=10)

## Trace plots of the parameters
par(mai=rep(0.5,4))
plot(spatial.logist$p.beta.theta.samples)


## Here we take the burn-in to be 90% of the iterations
burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

## Here we calculate the summary statistics for the samples of the
## covariance parameters after the burn-in
print(summary(window(spatial.logist$p.beta.theta.samples, start=burn.in)))

## Here we derive the probability of having malaria.
## For this we take the beta samples, the samples for the spatial random effects and
## we apply the formula to compute the predicted probability in a logistic regression
beta.hat <- spatial.logist$p.beta.theta.samples[sub.samps,1:3]
eta.hat <- spatial.logist$p.w.samples[,sub.samps]
p.hat <- matrix(0,dim(un.coord)[1],dim(beta.hat)[1])
for(k in 1:dim(beta.hat)[1]){
    p.hat[,k] <- exp(beta.hat[k,1]+beta.hat[k,2]*net.use.coord.gambia+beta.hat[k,3]*green.coord.gambia+eta.hat[,k])/(1+exp(beta.hat[k,1]+beta.hat[k,2]*net.use.coord.gambia+beta.hat[k,3]*green.coord.gambia+eta.hat[,k]))
}

## This is the estimated probability of malaria (at each location, we are taking the median
## of the predicted probability)
p.hat.median <- apply(p.hat,1,median)


## Here we estimate the spatial random effects and we compute the median, and the extremes of a 95% confidence interval
eta.post.median <- apply(eta.hat,1,median)
eta.post.low.bd <- apply(eta.hat,1,quantile,0.025)
eta.post.upp.bd <- apply(eta.hat,1,quantile,0.975)

## This is to make a plot of the estimated spatial random effects
surf.eta <- mba.surf(cbind(un.coord,eta.post.median),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image.plot(surf.eta, main="Interpolated posterior median \n of spatial random effects")
contour(surf.eta, add=TRUE)
points(un.coord[,1], un.coord[,2],pch=19,col="black")

## This is to make a plot of the lower bound of the 95% CI the estimated spatial random effects
surf.eta.low <- mba.surf(cbind(un.coord,eta.post.low.bd),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image.plot(surf.eta.low, main="Interpolated lower bound of  \n 95% CI for spatial random effects",zlim=c(-3,4))
contour(surf.eta.low, add=TRUE)
points(un.coord[,1], un.coord[,2],pch=19,col="black")

surf.eta.upp <- mba.surf(cbind(un.coord,eta.post.upp.bd),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image.plot(surf.eta.upp, main="Interpolated upper bound of \n of 95% CI for spatial random effects",zlim=c(-3,4))
contour(surf.eta.upp, add=TRUE)
points(un.coord[,1], un.coord[,2],pch=19,col="black")


## This is to make a plot of the estimated probability of malaria
surf.p <- mba.surf(cbind(un.coord,p.hat.median),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image.plot(surf.p, main="Interpolated posterior median of \n probability of malaria",col=heat.colors(100)[90:1])
contour(surf.p, add=TRUE)
points(un.coord[,1], un.coord[,2],pch=19,col="black")




