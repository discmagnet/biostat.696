# Homework 03 R Code
setwd("~/WORKING_DIRECTORIES/biostat.696/homework03")

# Import dataset
election <- read_csv("election_data_2016.csv")

# Load packages
library(maps)
library(maptools)
library(spdep)
library(classInt)
library(RColorBrewer)
library(SpatialEpi)

# Select just the continental United States
# Removes Alaska, District of Columbia, and Hawaii
election <- election[-c(2,9,12),]
US <- map("state",fill=TRUE,plot=FALSE)
US.names <- US$names
US.IDs <- sapply(strsplit(US.names,":"),function(x) x[1])
US_poly_sp <- map2SpatialPolygons(US,IDs=US.IDs,proj4string=CRS("+proj=longlat + datum=wgs84"))
index.columbia <- 8
US_no_columbia <- US_poly_sp[seq(1,49)[-8]]
us.nb <- poly2nb(US_no_columbia)
us.weights <- nb2WB(us.nb)
adj.us <- us.weights$adj
weights.us <- us.weights$weights
num.us <- us.weights$num
us.listw <- nb2listw(us.nb)
# Create choropleth map of percentage of Trump voters
plotvar <- election$pct.Trump
nclr <- 5
plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="quantile")
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Percentage of Trump Voters")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[31.9%, 41.16%)","[41.16%, 47.52%)","[47.52%, 52.98%)","[52.98%, 59.66%)","[59.66%, 70.1%]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
# Test for spatial correlation by looking at Moran's I statistic
moran.trump <- moran.test(election$pct.Trump,listw=us.listw,randomisation=FALSE)
moran.trump
# Select variables
election$income_c <- election$`median household income 2016`-57059
election$unemploy <- election$`unemployment rate`
election$white <- election$`pct white`
election$college <- election$pct_bachelor_degree
# Test for spatial correlation by looking at Moran's I of residuals of regression
lm.morantest(lm(election$pct.Trump~election$income_c+election$unemploy+election$white+election$college),
             listw=us.listw)
# Choropleth map of Median Household Income in 2016
plotvar <- election$`median household income 2016`
nclr <- 5
plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="quantile")
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Median Household Income in 2016")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[$41,754, $49,711)","[$49,711, $53,462)","[$53,462, $57,048)","[$57,048, $65,651)","[$65,651, $78,945]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
# Choropleth map of Unemployment Rate
plotvar <- election$`unemployment rate`
nclr <- 5
plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="quantile")
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Unemployment Rate")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[2.8%, 5.38%)","[5.38%, 6.68%)","[6.68%, 7.5%)","[7.75%, 8.26%)","[8.26%, 9.6%]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
# Choropleth map of Percentage who are White
plotvar <- election$`pct white`
nclr <- 5
plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="quantile")
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Percentage who are White")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[57.2%, 68.7%)","[68.7%, 77.26%)","[77.26%, 82.8%)","[82.8%, 87.8%)","[87.8%, 94.8%]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
# Choropleth map of Percentage who have a Bachelor's Degree
plotvar <- election$pct_bachelor_degree
nclr <- 5
plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="quantile")
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Percentage who have a Bachelor's Degree")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[11.9%, 16.18%)","[16.18%, 17.88%)","[17.88%, 19.44%)","[19.44%, 20.88%)","[20.88%, 24.4%]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
# Fit an improper CAR model
library(CARBayes)
formula <- election$pct.Trump~election$income_c+election$unemploy+election$white+election$college
N <- length(election$pct.Trump)
rep.us <- rep(1:N,num.us)
W <- matrix(0,N,N)
for(i in 1:N){
  W[i,adj.us[rep.us==i]] <- rep(1,num.us[i])
}
model.car <- S.CARleroux(formula=formula, W=W, family="gaussian",
                         rho=1,burnin=20000, n.sample=3000000,thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,prior.nu2=NULL,prior.tau2=NULL, verbose=TRUE)
model.car
samples.eta <- model.car$samples$phi
post.median.eta <- as.numeric(apply(samples.eta,2,median))
plotvar <- post.median.eta
nclr <- 4
plotclr <- brewer.pal(nclr,"RdBu")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(plotvar),summary(plotvar)[2],0,summary(plotvar)[5],max(plotvar)),length=nclr+1)
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Posterior Median of Spatial Random Effects")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[-17.42, -3.77)","[-3.77, 0)","[0, 3.19)","[3.19, 9.41]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
low.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.025))
plotvar <- low.bd95ci.eta
nclr <- 4
plotclr <- brewer.pal(nclr,"Reds")[4:1]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(plotvar),summary(plotvar)[2],summary(plotvar)[3],summary(plotvar)[5],max(plotvar)),length=nclr+1)
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Posterior Lower Bound of Spatial Random Effects")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[-22.62, -7.76)","[-7.76, -2.77)","[-2.77, -0.21)","[-0.21, -0.10]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
up.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.975))
plotvar <- up.bd95ci.eta
nclr <- 4
plotclr <- brewer.pal(nclr,"Blues")[1:4]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(plotvar),summary(plotvar)[2],summary(plotvar)[3],summary(plotvar)[5],max(plotvar)),length=nclr+1)
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Posterior Upper Bound of Spatial Random Effects")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[0.13, 0.22)","[0.22, 3.04)","[3.04, 8.66)","[8.66, 13.82]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
# Fit a proper CAR model
us.w.proper.car <- nb2listw(us.nb,style="B")
model.proper.car <- spautolm(formula=formula,listw=us.w.proper.car,family="CAR")
summary(model.proper.car)
# Fit a SAR model
model.sar <- spautolm(formula,listw=us.listw,family="SAR")
summary(model.sar)
# Disease mapping
election$voterthousand <- ceiling(election$Trump/1000)
formula2 <- election$voterthousand ~ election$income_c+election$unemploy+election$white+
  election$college+offset(log(election$No_voters))
model.car.poisson <- S.CARleroux(formula=formula2, W=W, family="poisson",
                                 rho=1,burnin=20000, n.sample=2500000,thin=1,
                                 prior.mean.beta=NULL, prior.var.beta=NULL,
                                 prior.nu2=NULL,prior.tau2=NULL, verbose=TRUE)
model.car.poisson$samples$beta
samples.eta <- model.car.poisson$samples$phi
fitted.val <- matrix(NA,dim(model.car.poisson$samples$beta)[1],48)
for(j in 1:dim(model.car.poisson$samples$beta)[1]){
  fitted.val[j,] <- model.car.poisson$samples$beta[j,1] +
    model.car.poisson$samples$beta[j,2]*election$income_c +
    model.car.poisson$samples$beta[j,3]*election$unemploy +
    model.car.poisson$samples$beta[j,4]*election$white +
    model.car.poisson$samples$beta[j,5]*election$college + 
    samples.eta[j,]
  print(j)
}
samples.eta <- model.car.poisson$samples$phi
post.median.eta <- as.numeric(apply(samples.eta,2,median))
plotvar <- post.median.eta
nclr <- 4
plotclr <- brewer.pal(nclr,"RdBu")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(plotvar),summary(plotvar)[2],0,summary(plotvar)[5],max(plotvar)),length=nclr+1)
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Posterior Median of Spatial Random Effects")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[-0.34, -0.10)","[-0.10, 0)","[0, 0.10)","[0.10, 0.41]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
low.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.025))
plotvar <- low.bd95ci.eta
nclr <- 4
plotclr <- brewer.pal(nclr,"Reds")[4:1]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(plotvar),summary(plotvar)[2],summary(plotvar)[3],summary(plotvar)[5],max(plotvar)),length=nclr+1)
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Posterior Lower Bound of Spatial Random Effects")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[-0.53, -0.21)","[-0.21, -0.10)","[-0.10, -0.03)","[-0.03, 0.28]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
up.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.975))
plotvar <- up.bd95ci.eta
nclr <- 4
plotclr <- brewer.pal(nclr,"Blues")[1:4]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(plotvar),summary(plotvar)[2],summary(plotvar)[3],summary(plotvar)[5],max(plotvar)),length=nclr+1)
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Posterior Upper Bound of Spatial Random Effects")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("[-0.17, 0.01)","[0.01, 0.11)","[0.11, 0.24)","[0.24, 0.53]")
legend(x=-75,y=40,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
plotvar <- -1*(up.bd95ci.eta < 0) + 1*(low.bd95ci.eta > 0)
nclr <- 3
plotclr <- brewer.pal(nclr,"RdBu")[3:1]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(-1.5,-0.5,0.5,1.5),length=nclr+1)
class
colcode <- findColours(class,plotclr)
plot(US_no_columbia,border="black",axes=TRUE,xlim=c(-125,-55))
title(xlab="Longitude",ylab="Latitude",main="Significance of Spatial Random Effects")
plot(US_no_columbia,col=colcode,add=T)
leg.txt<-c("Significantly Negative","Insignificant","Significantly Positive")
legend(x=-78,y=35,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")