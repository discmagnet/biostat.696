


########################################################################################
####
####     Bayesian hierarchical smoothing of the percent burglaries in Michigan counties in 2006.
####     The dataset contains information on the unemployment rate,
####     the median housing income, and the high school graduation rate
####     in Michigan counties in year 2006
####
########################################################################################


library(maps)
library(spdep)
library(classInt)
library(RColorBrewer)
library(SpatialEpi)
library(maptools)
library(nlme)
library(CARBayes)

## Reading the data in R
mi.crime <- read.csv("Crime_data_Michigan_2006.csv", sep=",",header=TRUE)

## These are the variables in the dataset
burglary <- mi.crime$Burglary
total.crime <- mi.crime$Total_no_offenses
county <- mi.crime$County
unempt <- mi.crime$Unemployment_rate
high.school <- mi.crime$High.school.graduate.or.higher_25to34
income <- mi.crime$Median_household_income
pct.burglary <- burglary/total.crime

## This is to get the neighboring structure for the counties in Michigan
mi.county<-map("county","michigan",fill=T,plot=F)
mi.IDs <- sapply(strsplit(mi.county$names,":"),function(x) x[1])
mi.poly.county <- map2SpatialPolygons(mi.county,IDs=mi.IDs,proj4string=CRS("+proj=longlat + datum=wgs84"))

mi.nb <- poly2nb(mi.poly.county)
mi.nb

## Here we get the adjacency weights, together with the information on the number of weights for each county
mi.weights <- nb2WB(mi.nb)
adj.mi <- mi.weights$adj
weights.mi <- mi.weights$weights
num.mi <- mi.weights$num
mi.listw <- nb2listw(mi.nb)


### This is some code to create the adjacency matrix W
N <- length(pct.burglary)
rep.mi <- rep(1:N,num.mi)
W <- matrix(0,N,N)
for(i in 1:N){
    W[i,adj.mi[rep.mi==i]] <- rep(1,num.mi[i])
}


## Making plots of the response variable
plotvar <- pct.burglary*100
nclr <- 5

range(plotvar)
#[1]  1.73913 12.87313

plotclr <- brewer.pal(nclr,"Greys")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(plotvar,c(0,0.2,0.4,0.6,0.8,1.0)))

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Percentage of burglaries \n Michigan, 2006")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[1.74%;4.63%)","[4.63%;5.41%)","[5.41%;6.42%)","[6.42%;8.18%)","[8.18%;12.87%]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


## Making plots of the covariates
# Unemployment rate
plotvar <- unempt
nclr <- 5

range(plotvar)
#[1]  4.8 13.3

plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(plotvar,c(0,0.2,0.4,0.6,0.8,1.0)))

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Umeployment rate \n Michigan, 2006")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[4.8;6.38)","[6.38;7.6)","[7.6;8.12)","[8.12;8.86)","[8.86;13.3]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


# High school
plotvar <- high.school
nclr <- 5

range(plotvar)
#[1] 64.4 96.5

plotclr <- brewer.pal(nclr,"Blues")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(plotvar,c(0,0.2,0.4,0.6,0.8,1.0)))
class
#[64.4,86.88) [86.88,88.78)  [88.78,90.5)  [90.5,93.04)  [93.04,96.5]

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="High school graduation rate \n Michigan, 2006")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[64.4;86.88)","[86.88;88.78)","[88.78;90.5)","[90.5;93.04)","[93.04;96.5]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")



# Income
plotvar <- income/1000
nclr <- 5

range(plotvar)
#[1]  3.447 72.129

plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(plotvar,c(0,0.2,0.4,0.6,0.8,1.0)))
class
#[3.447,37.1082) [37.1082,40.178)  [40.178,43.271) [43.271,48.9536) [48.9536,72.129]

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Median annual housing income \n in thousands of dollars")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[3.45;37.11)","[37.11,40.18)","[40.18,43.27)","[43.27,48.95)","[48.95,72.13]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")




## Here we are computing Moran's I to see if the data is spatially correlated
moran.test(pct.burglary, mi.listw, rank=TRUE)


### Fitting a linear model to the percent burglary with median housing income (in 1000 of $)
### unemployment rate, and high school graduation rate
income.thous <- income/1000

lm.burglary <- lm(pct.burglary~income.thous+high.school+unempt)
summary(lm.burglary)

Coefficients:
Estimate Std. Error t value Pr(>|t|)
(Intercept)   0.0795126  0.0623766   1.275   0.2061
income.thous -0.0001215  0.0003643  -0.333   0.7398
high.school  -0.0004308  0.0005853  -0.736   0.4638
unempt        0.0035835  0.0018861   1.900   0.0611 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.02257 on 79 degrees of freedom
Multiple R-squared:  0.1266,    Adjusted R-squared:  0.09342



## Computing Moran's I on the residuals to determine if there is residual
## spatial variation
res.burglary <- as.numeric(summary(lm.burglary)$residuals)
moran.test(res.burglary, mi.listw, rank=TRUE)


######################################################################
######################################################################
###
###
### This code fits a Bayesian smoothing model to the percent burglaries
### using an improper CAR prior for the spatial random effects.
### We are using the CARBayes package, the S.CARleroux function and all the default choices
### for that function
###
######################################################################
######################################################################

formula <- pct.burglary~income.thous+high.school+unempt
model.car <- S.CARleroux(formula=formula, W=W, family="gaussian",
rho=1,burnin=10000, n.sample=40000,thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,prior.nu2=NULL,prior.tau2=NULL, verbose=TRUE)


model.car$summary.results
Median    2.5%  97.5% n.sample % accept n.effective Geweke.diag
(Intercept)   0.0858 -0.0942 0.2658    30000      100     11217.7         0.5
income.thous -0.0007 -0.0018 0.0004    30000      100      7253.8         0.0
high.school  -0.0001 -0.0018 0.0015    30000      100     12971.6        -0.8
unempt        0.0026 -0.0029 0.0082    30000      100      9592.9         0.6
nu2           0.0007  0.0005 0.0010    30000      100     14607.4         0.0
tau2          0.0012  0.0008 0.0021    30000      100      6356.7         0.2
rho           1.0000  1.0000 1.0000       NA       NA          NA          NA


#####################################################################
#####################################################################
###
###       SPATIAL PLOTS
###
### Here we are going to make spatial maps of the posterior
### medians and posterior standard deviations
### of the spatial effects, and of the smoothed fitted values
###
###
#####################################################################
#####################################################################

names(model.car)
[1] "summary.results"     "samples"             "fitted.values"       "residuals"
"modelfit"            "accept"              "localised.structure" "formula"
[9] "model"               "X"


### Posterior median of spatial effects
samples.eta <- model.car$samples$phi
post.median.eta <- as.numeric(apply(samples.eta,2,median))
plotvar <- post.median.eta
nclr <- 5

plotclr <- brewer.pal(nclr,"RdBu")[5:1]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(-0.045,-0.02,-0.004,0,0.0,0.006,0.017),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Spatial random effects \n CAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[-0.04;-0.02)","[-0.02;-0.004)","[-0.004;0.0)","[0.0;0.006)","[0.006;0.017]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


### Posterior SD of spatial effects
post.sd.eta <- as.numeric(apply(samples.eta,2,sd))
plotvar <- post.sd.eta
nclr <- 5

plotclr <- brewer.pal(nclr,"YlGnBu")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(0.012,0.013,0.014,0.016,0.017,0.023),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="SD of spatial random effects \n CAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[0.012; 0.013)","[0.013; 0.014)","[0.014; 0.016)","[0.016; 0.017)","[0.017; 0.023]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")



### Lower bd 95% CI of spatial effects
low.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.025))
plotvar <- low.bd95ci.eta
nclr <- 5

plotclr <- brewer.pal(nclr,"Blues")[5:1]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(-0.082,-0.049,-0.032,-0.027,-0.023,-0.011),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Lower bound of 95% CI for spatial random effects \n CAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[-0.082; -0.049)","[-0.049; -0.032)","[-0.032; -0.027)","[-0.027; -0.023)","[-0.023; -0.011]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


### Upper bd 95% CI of spatial effects
upp.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.975))
plotvar <- upp.bd95ci.eta
nclr <- 5


plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(-0.003,0.0,0.025,0.027,0.035,0.057),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Upper bound of 95% CI for spatial random effects \n CAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[-0.003; 0.0)","[0.017; 0.025)","[0.025; 0.027)","[0.027; 0.035)","[0.035; 0.057]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")




### Posterior median of fitted values
fitted.val <- matrix(NA,dim(model.car$samples$beta)[1],83)
for(j in 1:dim(model.car$samples$beta)[1]){
    fitted.val[j,] <- model.car$samples$beta[j,1] +
    model.car$samples$beta[j,2]*income.thous +
    model.car$samples$beta[j,3]*high.school +
    model.car$samples$beta[j,4]*unempt + samples.eta[j,]
    print(j)
}


post.median.fitted <- as.numeric(apply(fitted.val,2,median))
plotvar <- post.median.fitted*100
nclr <- 5

plotclr <- brewer.pal(nclr,"Greys")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(pct.burglary*100,c(0.0,0.2,0.4,0.6,0.8,1.0)),length=nclr+1)


colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Fitted values \n CAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[1.74%;4.63%)","[4.63%;5.41%)","[5.41%;6.42%)","[6.42%;8.18%)","[8.18%;12.87%]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")




######################################################
###
### Here we are fitting a SAR model to the data
###
###
######################################################


model.sar <- spautolm(pct.burglary~income.thous+high.school+unempt,listw=mi.listw,family="SAR")
summary(model.sar)


Coefficients:
Estimate  Std. Error z value Pr(>|z|)
(Intercept)   0.06907423  0.05903469  1.1701   0.2420
income.thous -0.00047499  0.00036906 -1.2870   0.1981
high.school  -0.00013117  0.00054232 -0.2419   0.8089
unempt        0.00339504  0.00179022  1.8964   0.0579

Lambda: 0.36475 LR test value: 4.4144 p-value: 0.035637
Numerical Hessian standard error of lambda: 0.15755

Log likelihood: 201.1588
ML residual variance (sigma squared): 0.00044565, (sigma: 0.02111)
Number of observations: 83
Number of parameters estimated: 6
AIC: -390.32



### Fitted values according to a SAR model
fitted.val.sar <- as.numeric(model.sar$fit$fitted.values)
plotvar <- fitted.val.sar*100
nclr <- 5

plotclr <- brewer.pal(nclr,"Greys")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(pct.burglary*100,c(0.0,0.2,0.4,0.6,0.8,1.0)),length=nclr+1)

colcode <- findColours(class,plotclr)


plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Fitted values \n SAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[1.74%;4.63%)","[4.63%;5.41%)","[5.41%;6.42%)","[6.42%;8.18%)","[8.18%;12.87%]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")



######################################################
###
### Here we are fitting a proper CAR model to the data
###
###
######################################################

mi.w.proper.car <- nb2listw(mi.nb,style="B")

model.proper.car <- spautolm(pct.burglary~income.thous+high.school+unempt,listw=mi.w.proper.car,family="CAR")
summary(model.proper.car)

Coefficients:
Estimate  Std. Error z value Pr(>|z|)
(Intercept)   0.06521504  0.06028501  1.0818   0.2794
income.thous -0.00036177  0.00037273 -0.9706   0.3318
high.school  -0.00015984  0.00055334 -0.2889   0.7727
unempt        0.00332211  0.00181305  1.8323   0.0669

Lambda: 0.10513 LR test value: 2.289 p-value: 0.13029
Numerical Hessian standard error of lambda: 0.057572

Log likelihood: 200.0961
ML residual variance (sigma squared): 0.00045485, (sigma: 0.021327)
Number of observations: 83
Number of parameters estimated: 6
AIC: -388.19


names(model.proper.car)
[1] "fit"         "lambda"      "LL"          "LL0"         "call"        "parameters"  "aliased"     "method"      "family"
[10] "zero.policy" "weights"     "interval"    "trs"         "timings"     "LLNullLlm"   "fdHess"      "lambda.se"   "X"
[19] "Y"

### Fitted values according to a proper CAR model
fitted.val.proper.car <- as.numeric(model.proper.car$fit$fitted.values)
plotvar <- fitted.val.proper.car*100
nclr <- 5

plotclr <- brewer.pal(nclr,"Greys")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(pct.burglary*100,c(0.0,0.2,0.4,0.6,0.8,1.0)),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Fitted values \n proper CAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[1.74%;4.63%)","[4.63%;5.41%)","[5.41%;6.42%)","[6.42%;8.18%)","[8.18%;12.87%]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")




