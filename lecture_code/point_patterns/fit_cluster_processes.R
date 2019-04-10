
######  This code shows how to fit a inhomogenous Poisson process and a cluster point
######  process to the data



library(splancs)
library(spatstat)
library(maptools)
library(spdep)
library(spatialkernel)


## Data on Beilschmieda trees

data(bei)
bei
Planar point pattern: 3604 points
window: rectangle = [0, 1000] x [0, 500] metres

plot(bei,main="Beilschmieda tree")


## Fitting a homogenous Poisson process to this data
fit <- ppm(bei~1)
fit
Stationary Poisson process
Intensity: 0.007208
Estimate       S.E.   CI95.lo   CI95.hi Ztest      Zval
log(lambda) -4.932564 0.01665742 -4.965212 -4.899916   *** -296.1182


### Fitting an inhomogenous Poisson process using the information in the dataset bei.extra
### which contains information on the terrain elevation and slope
bei.extra
List of pixel images

elev:
real-valued pixel image
101 x 201 pixel array (ny, nx)
enclosing rectangle: [-2.5, 1002.5] x [-2.5, 502.5] metres

grad:
real-valued pixel image
101 x 201 pixel array (ny, nx)
enclosing rectangle: [-2.5, 1002.5] x [-2.5, 502.5] metres

par(mfrow=c(1,1))
plot(bei.extra$elev, main="Elevation",col=terrain.colors(100))
par(mfrow=c(1,1))
plot(bei.extra$grad, main="Terrain slope",col=terrain.colors(100))


## Modeling the log intensity function as a function of the terrain slope
fit <- ppm(bei~grad, data=bei.extra)
fit

Nonstationary Poisson process

Log intensity:  ~grad

Fitted trend coefficients:
(Intercept)        grad
-5.390553    5.022021

Estimate       S.E.   CI95.lo   CI95.hi Ztest       Zval
(Intercept) -5.390553 0.03001716 -5.449385 -5.331720   *** -179.58236
grad         5.022021 0.24540264  4.541041  5.503002   ***   20.46441


### Plotting how the intensity changes as a function of terrain slope
plot(effectfun(fit,"grad",se.fit=TRUE),main="Estimated intensity with 95% CI",xlab="Terrain slope")


## Modeling the log intensity function as a quadratic function of terrain slope
fit <- ppm(bei~grad+I(grad^2), data=bei.extra)
fit
Nonstationary Poisson process

Log intensity:  ~grad + I(grad^2)

Fitted trend coefficients:
(Intercept)        grad   I(grad^2)
-5.987136   18.744888  -55.602328

Estimate       S.E.    CI95.lo    CI95.hi Ztest       Zval
(Intercept)  -5.987136 0.05573988  -6.096384  -5.877888   *** -107.41207
grad         18.744888 1.06264319  16.662146  20.827630   ***   17.63987
I(grad^2)   -55.602328 4.27599072 -63.983116 -47.221540   ***  -13.00338


### Plotting how the intensity changes as a function of terrain slope
plot(effectfun(fit,"grad",se.fit=TRUE),main="Estimated intensity with 95% CI \n Quadratic polynomial",xlab="Terrain slope")


### Using a factor covariate
data(gorillas)
gorillas
Marked planar point pattern: 647 points
Mark variables: group, season, date
window: polygonal boundary
enclosing rectangle: [580457.9, 585934] x [674172.8, 678739.2] metres

## Rescaling the window to have unit of measure in km rather than meters
gor <- rescale(gorillas, 1000, unitname="km")
## We only care about the unmarked data (e.g. the locations of the gorilla nests)
gor <- unmark(gor)

plot(gor,main="Gorilla nests",pch=19)


## The dataset gorillas.extra has extra information, including the type of vegetation or cover type
gorillas.extra
List of pixel images

aspect:
factor-valued pixel image
factor levels:
[1] "N"  "NE" "E"  "SE" "S"  "SW" "W"  "NW"
149 x 181 pixel array (ny, nx)
enclosing rectangle: [580440, 586000] x [674160, 678730] metres

elevation:
integer-valued pixel image
149 x 181 pixel array (ny, nx)
enclosing rectangle: [580440, 586000] x [674160, 678730] metres

### Here we rescale also the dataset with the extra information using the same unit
### of scale that we used on the data with the gorillas nests
gex <- lapply(gorillas.extra,rescale,s=1000,unitname="km")

names(gex)
[1] "aspect"     "elevation"  "heat"       "slopeangle" "slopetype"  "vegetation" "waterdist"

### We only care the vegetation factor for the moment
levels(gex$vegetation)
[1] "Disturbed"  "Colonising" "Grassland"  "Primary"    "Secondary"  "Transition"

## Since the names are too long, we use this function to shorten each name to
## only the first 4 letters
shorten <- function(x){substr(x,1,4)}
names(gex) <- shorten(names(gex))
names(gex)
names(gex)[4:5] <- c("sang","styp")

isfactor <- !sapply(lapply(gex,levels),is.null)
for(i in which(isfactor)){
    levels(gex[[i]]) <- shorten(levels(gex[[i]]))
}

names(gex)
levels(gex$vege)
## Here we plot the vegetation types, since there are only 6 types we only use 6 colors
plot(gex$vege,col=terrain.colors(6),main="Type of vegetation or cover")


### Fitting an inhomogenous Poisson process with vegetation as covariate
fit.vege1 <- ppm(gor~vege, data=gex)
fit.vege1
Nonstationary Poisson process

Log intensity:  ~vege

Fitted trend coefficients:
(Intercept)    vegeColo    vegeGras    vegePrim    vegeSeco    vegeTran
2.3238438   2.0816609  -0.7720732   2.1147512   1.1341176   1.6151353

            Estimate      S.E.    CI95.lo    CI95.hi Ztest      Zval
(Intercept)  2.3238438 0.1059998  2.1160881  2.5315996   *** 21.923099
vegeColo     2.0816609 0.5870002  0.9311615  3.2321602   ***  3.546269
vegeGras    -0.7720732 0.2474590 -1.2570840 -0.2870625    ** -3.120005
vegePrim     2.1147512 0.1151001  1.8891592  2.3403432   *** 18.373152
vegeSeco     1.1341176 0.2426005  0.6586294  1.6096058   ***  4.674836
vegeTran     1.6151353 0.2646875  1.0963573  2.1339133   ***  6.102046


fit.vege2 <- ppm(gor~vege-1,data=gex)
fit.vege2
Nonstationary Poisson process

Log intensity:  ~vege - 1

Fitted trend coefficients:
vegeDist vegeColo vegeGras vegePrim vegeSeco vegeTran
2.323844 4.405505 1.551771 4.438595 3.457961 3.938979

        Estimate       S.E.  CI95.lo  CI95.hi Ztest      Zval
vegeDist 2.323844 0.10599979 2.116088 2.531600   *** 21.923099
vegeColo 4.405505 0.57735027 3.273919 5.537090   ***  7.630558
vegeGras 1.551771 0.22360680 1.113509 1.990032   ***  6.939729
vegePrim 4.438595 0.04485613 4.350679 4.526511   *** 98.951803
vegeSeco 3.457961 0.21821789 3.030262 3.885661   *** 15.846370
vegeTran 3.938979 0.24253563 3.463618 4.414340   *** 16.240827


### Since the intensity is constant within each category, we
### can obtain what is the intensity within each type of vegetation type
exp(coef(fit.vege2))
vegeDist vegeColo vegeGras vegePrim vegeSeco vegeTran
10.21486 81.90047  4.71982 84.65592 31.75218 51.36614


### Fitting a model for the intensity function that depends only on the coordinates
fit.coord <- ppm(gor~x+y)
fit.coord
Nonstationary Poisson process

Log intensity:  ~x + y

Fitted trend coefficients:
(Intercept)            x            y
-118.8650482   -0.3622621    0.4927706

Estimate        S.E.      CI95.lo     CI95.hi Ztest       Zval
(Intercept) -118.8650482 29.40243108 -176.4927542 -61.2373422   ***  -4.042695
x             -0.3622621  0.03040420   -0.4218533  -0.3026710   *** -11.914869
y              0.4927706  0.03736216    0.4195421   0.5659991   ***  13.189029


#### Modeling the intensity function using an offset
data(chorley)
chorley
Marked planar point pattern: 1036 points
Multitype, with levels = larynx, lung
window: polygonal boundary
enclosing rectangle: [343.45, 366.45] x [410.41, 431.79] km

## We are splitting the lung and larynx cancer
lung <- split(chorley)$lung
larynx <- split(chorley)$larynx


plot(lung$x,lung$y,main="Lung cancer cases",col="red",pch=19,xlab="Easting",ylab="Northing")
plot(larynx$x,larynx$y,main="Larynx cancer cases",col="blue",pch=19,xlab="Easting",ylab="Northing")


## We would like to have an idea of the population density.
## As we don't have information on this, we proceed by taking the lung cases and
## smoothing it out, assuming that this can give us a good idea of the spatially-varying population density.
## This does kernel density smoothing using a bandwidth of 0.15 with pixels that are 0.1km wide.
## Additionally since the density has to be positive, negative or zero values are replaced with small
## positive numbers
smo <- density(lung,sigma=0.15,eps=0.1,positive=TRUE)
smo <- eval.im(pmax(smo,1e-10))

plot(smo,main="Smooth intensity for lung cancer cases \n or approx. population density")

## Now we use the smooth population density as offset and model the
## the intensity for the larynx cancer cases as
## lambda(s)= beta*pop(s)
fit.pop <- ppm(larynx~offset(log(smo)))
fit.pop
Nonstationary Poisson process

Log intensity:  ~offset(log(smo))

Fitted trend coefficient:  (Intercept) = -2.938568

Estimate      S.E.   CI95.lo   CI95.hi Ztest      Zval
(Intercept) -2.938568 0.1313064 -3.195924 -2.681212   *** -22.37947

#######

# Fitting a Cox process model to the data via minimum contrast.
# The intensity is modeled as a log Gaussian process with exponential correlation function.
data(redwood)

# The function lgcp.estK fit a log Gaussian Cox Process to the data
# by computing the K function for the data and comparing it to the theoretical K function
# for a log Gaussian Cox process.
# The function only works for intensity processes that are log Gaussian processes, stationary and
# isotropic with an exponential
# covariance function.
# In order to compute the log Gaussian Cox Process, it is necessary to specify the parameters of
# the exponential covariance function:
# the (partial) sill sigma2, and the range parameter (phi), here called alpha
# Then we need to specify the parameters for the least square procedure
# (here called method of minimum contrasts), that is the power c, here called
# q, a second power term p, that for least squares should be set equal to 2,
# and the bounds for the integral, rmin and rmax.
# Default values are: c (that is, q) = 1/4, p=2 and rmin and rmax the ranges for r used by the Kest function

u.logcp <- lgcp.estK(redwood, c(sigma2=0.1, alpha=1))
u.logcp

Minimum contrast fit (object of class “minconfit”)
Model: log-Gaussian Cox process
Fitted by matching theoretical K function to Kest(redwood)
Parameters fitted by minimum contrast ($par):
sigma2     alpha
1.0485493 0.0997963
Derived parameters of log-Gaussian Cox process ($modelpar):
sigma2     alpha        mu
1.0485493 0.0997963 3.6028597
Converged successfully after 145 iterations.
Domain of integration: [ 0 , 0.25 ]
Exponents: p= 2, q= 0.25

plot(u.logcp,main="Fitted K function and theoretical K function \n Redwood data, log Gaussian-Cox process")

lty col    key       label                                      meaning
fit      1   1    fit   K[fit](r)                 minimum contrast fit of LGCP
iso      2   2    iso   K[iso](r) Ripley isotropic correction estimate of K(r)
trans    3   3  trans K[trans](r)       translation-corrected estimate of K(r)
border   4   4 border  K[bord](r)            border-corrected estimate of K(r)
theo     5   5   theo  K[pois](r)                     theoretical Poisson K(r)



# Other models for which least square fitting is provided is a Poisson cluster process (matclust.estK)

# Here we fit a Matern cluster process to the redwood data using the matclust.estK function
# The function requires as argument the starting values for the parameters of the Matern cluster process:
# the intensity of the parent process
# and the radius for the dispersal function. The parameter of the Poisson distribution for the offspring
# is derived from the estimated intensity of
# the entire process.
# Other arguments for the matclust.estK function include an optional argument, the intensity of the
# process, and the parameters q, p, rmin and rmax
# used in the least squares method.
u.matclust <- matclust.estK(redwood, c(kappa=10, R=0.1))
u.matclust

Minimum contrast fit (object of class “minconfit”)
Model: Matern Cluster process
Fitted by matching theoretical K function to Kest(redwood)
Parameters fitted by minimum contrast ($par):
kappa           R
24.56179311  0.08653217
Derived parameters of Matern Cluster process ($modelpar):
kappa           R          mu
24.56179311  0.08653217  2.52424567
Converged successfully after 83 iterations.
Domain of integration: [ 0 , 0.25 ]
Exponents: p= 2, q= 0.25

plot(u.matclust,main="Fitted K function and theoretical K function \n Redwood data, Matern cluster process")









