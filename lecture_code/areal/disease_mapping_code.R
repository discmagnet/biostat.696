

###  This file shows how to use the package CARBayes to fit a
###  model to areal data with a random effect for the area
###  modeled as an improper CAR random effect

###  We are going to use the data in the SpatialEpi package
###  that refers to Scotland and the number of lip cancer cases 
###  in the 56 counties of Scotland

###  We are modeling the number of observed cases for each county, O_i, 
###  as a Poisson random variable with parameter mu_i defined as follows.
###  The log risk for county i can be thought as the sum of a baseline log risk
###  alpha0 and a relative log risk for county i which is due to the effect of 
###  a covariate and a random effect for county i



library(SpatialEpi)
library(RColorBrewer)
library(classInt)
library(spdep)
library(maps)
library(maptools)

data(scotland)
names(scotland)

scotland$data[1:3,]

scotland.data <- scotland$data

Y <- scotland.data$cases 
E <- scotland.data$expected
X <- scotland.data$AFF

N <- length(Y)

### In order to fit a model with improper CAR random effects we need to know the list of adjacent county to each county in Scotland.
### We can get an automatic list of neighbors using the spdep package in R and the function
### poly2nb function.
### In order to do this, however, the data needs to be in a particular format: it needs to be 
### an object of the type SpatialPolygon.
### We can transform the data into the SpatialPolygon format using the SpatialEpi package and the function
### polygon2spatial.polygon.

scotland.polygon <- scotland$polygon$polygon
scotland.nrepeats <- scotland$polygon$nrepeats
scotland.names <- scotland$data$county.names
scotland.spatial.polygon <- polygon2spatial_polygon(scotland.polygon,coordinate.system="+proj=utm",
scotland.names,scotland.nrepeats)

### This makes a map of Scotland with the various counties
plot(scotland.spatial.polygon,border="black")

### By doing scotland.spatial.polygon[#i], 
### you can access the information for the i-th county in scotland
### Here we are coloring the first county red 
plot(scotland.spatial.polygon[1],col="red",add=T)


### Now that we have done the transformation, we can use the spdep package and get 
### the list of neighbors
### for our data using the function poly2nb

scotland.nb <- poly2nb(scotland.spatial.polygon)
scotland.nb


### This makes a map of scotland with a link between counties that are adjacent
scotland.coord <- coordinates(scotland.spatial.polygon)
plot(scotland.spatial.polygon,border="black")
plot(scotland.nb,scotland.coord,col="red",pch=20,add=TRUE)


### Having generated the poly2nb object, using the nb2WB (which stands for 
### neighborhood objects to WinBUGS), we get
### three vectors with information on the list of neighbors, the weights (all equal to 1)
### and the number of neighbors for each county
scotland.weights <- nb2WB(scotland.nb)
adj <- scotland.weights$adj
weights <- scotland.weights$weights
num <- scotland.weights$num


### Looking at the num file, we see that we have 3 regions that do not have 
### neighbors: these are the three
### islands. We can make two decisions at this point. Either we remove the 
### three islands from the analysis (the
### most pragmatic decision, based on the definition of neighborhoods) or we 
### include them in the analysis (but then we
### should use a different definition for neighborhoods based on maybe distance, 
### rather than contiguity)


#********************   Removing the three islands *****************************#

library(CARBayes)

### This is to identify the three regions in the dataset corresponding to 
### the three islands
index.islands <- which(num==0)

### To remove the islands from our spatial polygon data frame, we simply 
### remove the components
### of the list that correspond to the index for the three islands
scotland.noisland.spatial.polygon <- scotland.spatial.polygon[seq(1:N)[-index.islands]]

### Now we can apply the same commands we used earlier to get the list of neighbors, ### the weights
### and the number of neighbors that need to be passed to WinBUGS
plot(scotland.noisland.spatial.polygon,border="black")
plot(scotland.noisland.spatial.polygon[1],col="red",add=T)


### Here we get the list of neighbors for the data without islands using again 
### the function poly2nb

scotland.noisland.nb <- poly2nb(scotland.noisland.spatial.polygon)
scotland.noisland.nb

### This makes a map of scotland with a link between counties that are adjacent
scotland.noisland.coord <- coordinates(scotland.noisland.spatial.polygon)
plot(scotland.noisland.spatial.polygon,border="black")
plot(scotland.noisland.nb,scotland.noisland.coord,col="red",pch=20,add=TRUE)

### We then use the nb2WB function and we get
### the three vectors that we need to fit an improper CAR model in WinBUGS
scotland.noisland.weights <- nb2WB(scotland.noisland.nb)
adj.noisland <- scotland.noisland.weights$adj
weights.noisland <- scotland.noisland.weights$weights
num.noisland <- scotland.noisland.weights$num


### Since we have removed the three islands, we have to change also the data
### and remove from all the data the observations for the three islands

Y.noisland <- Y[-index.islands]
E.noisland <- E[-index.islands]
X.noisland <- X[-index.islands]

N.noisland <- length(Y.noisland)

### Now we fill the adjacency matrix W using the following trick 
W <- matrix(0,N.noisland,N.noisland)

rep.scotland.noisland <- rep(1:N.noisland,num.noisland)

for(i in 1:N.noisland){
	W[i,adj.noisland[rep.scotland.noisland==i]] <- rep(1,num.noisland[i])	
}


### The function S.CARbym fits spatial model to the Poisson data with two sets of random effects, the spatial
### and non spatial ones.
### Specifically, the model that will be fit is a Poisson model for the observed counts
### with mean (parameter) equal to E*relative risk where E represents the expected counts.
### In the S.CARbym function, we specify a formula that describes how the log of the
### mean of the Poisson distribution for the observed counts depends on covariates.
### Taking the log of E*relative risk, we have log(E)+log(relative risk).
### So the formula command expresses how the log(relative risk) depends on covariates and adds as offset
### the log(E), or the log of the expected counts.

### Here, since we assume that the log relative risk depends on
### the percentage of population
### involved in Agriculture, Fisheries and Forestries, we have:
formula <- Y.noisland~X.noisland+offset(log(E.noisland))

### Finally, the command that fits the Poisson model to the data is S.CARbym.
### The function takes as input
### the formula for the log relative risk that we just specified above,
### the adjacency matrix,
### the number of iterations for the MCMC algorithm (n.sample) and the number
### of iterations discarded for burnin (burnin).
model.scotland <- S.CARbym(formula=formula, family="poisson", W=W, burnin=5000, n.sample=10000)


### To see what the command poisson.iarCAR returns, just type
names(model.scotland)
[1] "summary.results"     "samples"             "fitted.values"
"residuals"           "modelfit"            "accept"
[7] "localised.structure" "formula"             "model"               "X"

names(model.scotland$samples)
[1] "beta"   "psi"    "tau2"   "sigma2" "fitted" "Y"


###
###
###       SPATIAL PLOTS
###
### Here we are going to make spatial maps of the posterior medians and posterior standard deviations
### of the spatial effects


### Raw data
plotvar <- Y.noisland
nclr <- 5

plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(plotvar,c(0,0.2,0.4,0.6,0.8,1.0)))
class

### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)


plot(scotland.noisland.spatial.polygon,border="black",axes=TRUE,xlim=c(-150,550))
title(xlab="Eastings (km)",ylab="Northings (km)",main="Observed number of cases")
plot(scotland.noisland.spatial.polygon,col=colcode,add=T)

leg.txt<-c("0-2","3-6","7-8","9-14",
"15-39")
legend(-150,1000,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")



### SMR or raw relative risks
plotvar <- Y.noisland/E.noisland
nclr <- 5

plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(plotvar,c(0,0.2,0.4,0.6,0.8,1.0)))
class

### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)


plot(scotland.noisland.spatial.polygon,border="black",axes=TRUE,xlim=c(-150,550))
title(xlab="Eastings (km)",ylab="Northings (km)",main="Standardized mortality ratio")
plot(scotland.noisland.spatial.polygon,col=colcode,add=T)

leg.txt<-c("[0; 0.37)","[0.37; 0.89)","[0.89; 1.23)","[1.23; 2.33)",
"[2.33; 6.43]")
legend(-150,1000,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

### AFF covariate
plotvar <- X
nclr <- 5

plotclr <- brewer.pal(nclr,"Blues")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(plotvar,c(0,0.2,0.4,0.6,0.8,1.0)))
class

### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)


plot(scotland.noisland.spatial.polygon,border="black",axes=TRUE,xlim=c(-150,550))
title(xlab="Eastings (km)",ylab="Northings (km)",main="Percentage of population involved in \n Agriculture, Fisheries and Forestries")
plot(scotland.noisland.spatial.polygon,col=colcode,add=T)

leg.txt<-c("[0%; 1%)","[1%; 7%)","[7%; 10%)","[10%; 16%)",
"[16%; 24%]")
legend(-150,1000,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


### Posterior median of random effects 
### The random.effects output from the S.CARbym function is matrix with N rows (where N is the number of areal units)
### and 5 columns: a column with the posterior mean, the posterior standard deviation, the posterior median, the lower 
### bound of a 95% credible interval and the upper bound of a 95% credible interval
post.median.psi <- as.numeric(apply(model.scotland$samples$psi,2,median))
plotvar <- post.median.psi
nclr <- 5

plotclr <- brewer.pal(nclr,"YlOrRd")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(plotvar,c(0,0.2,0.4,0.6,0.8,1.0)))
class

### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)

plot(scotland.noisland.spatial.polygon,border="black",axes=TRUE,xlim=c(-150,550))
title(xlab="Eastings (km)",ylab="Northings (km)",main="Estimated random effects")
plot(scotland.noisland.spatial.polygon,col=colcode,add=T)

leg.txt<-c("[-0.73,-0.52)","[-0.52,-0.29)","[-0.29,0.04)","[0.04,0.62)","[0.62,1.16]")
legend(-150,1000,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


### Estimated relative risks
### The relative risks are given by exp(X*beta+random effects) and are stored in the fitted values
 
###  This gets the component of the log relative risk that is given by X*beta.
###  The result of this is a matrix with as many rows as MCMC samples we have for the beta coefficients and as many columns
###  as number of areal units
RR.covariate <- model.scotland$samples$beta%*%t(matrix(cbind(rep(1,N.noisland),X.noisland),N.noisland,2))

###  This gets the component of the log relative risk that is given by the random effects
###  The result of this is a matrix with as many rows as MCMC samples we have for the spatial random effects (phi) and as many columns
###  as number of areal units

RR.random <- model.scotland$samples$psi
#+model.scotland$residuals[,3]
RR.all <- exp(RR.covariate+RR.random)

post.median.RR <- apply(RR.all,2,median)
plotvar <- post.median.RR
nclr <- 5

plotclr <- brewer.pal(nclr,"Blues")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(plotvar,c(0,0.2,0.4,0.6,0.8,1.0)))
class

### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)

plot(scotland.noisland.spatial.polygon,border="black",axes=TRUE,xlim=c(-150,550))
title(xlab="Eastings (km)",ylab="Northings (km)",main="Estimated relative risks")
plot(scotland.noisland.spatial.polygon,col=colcode,add=T)

leg.txt<-c("[0.35,0.51)","[0.51,0.83)","[0.83,1.18)","[1.18,1.96)","[1.96,4.54)")
legend(-150,1000,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")





