

######  This code shows how to estimate the intensity function
######  of a spatial point process non parametrically.


library(splancs)
library(spatstat)
library(maptools)
library(spdep)
library(spatialkernel)


data(japanesepines)
summary(japanesepines)


########################################################

## Here we estimate the intensity function by kernel smoothing.
## We use the density function in the spatstat package that only uses the Gaussian kernel
## and the bandwidth is then the standard deviation of the Gaussian kernel
int.japanese.01 <- density(japanesepines, 0.01)
plot(int.japanese.01,main="Kernel estimate of intensity function; h=0.01 \n Japanese black pines",xlab="X",ylab="Y")

int.japanese.05 <- density(japanesepines, 0.05)
plot(int.japanese.05,main="Kernel estimate of intensity function; h=0.05 \n Japanese black pines",xlab="X",ylab="Y")

int.japanese.1 <- density(japanesepines, 0.1)
plot(int.japanese.1,main="Kernel estimate of intensity function; h=0.1 \n Japanese black pines",xlab="X",ylab="Y")

int.japanese.5 <- density(japanesepines, 0.5)
plot(int.japanese.5,main="Kernel estimate of intensity function; h=0.5 \n Japanese black pines",xlab="X",ylab="Y")

## This gives the value of the intensity evaluated at points
density(japanesepines, 0.05, at="points")
density(japanesepines, 0.1, at="points")
density(japanesepines, 0.5, at="points")


################  Selecting the bandwidth

### To choose the correct bandwidth, we use the function mse2d in the splancs package
### In order to use that function we need to set the data in the "points" format
japanese_point <- matrix(NA, nrow=length(japanesepines$x),ncol=2)
japanese_point[,1] <- japanesepines$x
japanese_point[,2] <- japanesepines$y

## Here we plot the points
plot(japanese_point,main="Japanese black pines", xlab="X",ylab="Y",pch=20)

## Here we define a polygon that contains the events
japanese.poly <- as.points(c(min(japanese_point[,1]),max(japanese_point[,1]),max(japanese_point[,1]),min(japanese_point[,1])),
c(max(japanese_point[,2]),max(japanese_point[,2]),min(japanese_point[,2]),min(japanese_point[,2])))
plot(japanese.poly, type="l")


## This function computes the MSE for the estimated intensity function using different bandwidths
mse.japanese <- mse2d(japanese_point, poly=japanese.poly, nsmse=100, range=.5)

## Here we determine the optimal bandwidth as the one that minimizes the MSE
bw <- mse.japanese$h[which.min(mse.japanese$mse)]
bw

plot(x=mse.japanese$h, y=mse.japanese$mse, xlab="Bandwidth", ylab="MSE", type="l",main="MSE as a function of bandwidth h \n Japanese black pines")
i<-which.min(mse.japanese$mse)
points(mse.japanese$h[i], mse.japanese$mse[i])


## Here we plot the kernel estimate of the intensity function based on the selected bandwidth
int.japanese.bw <- density(japanesepines, bw)
plot(int.japanese.bw,main="Kernel estimate of intensity function; h=0.315 \n Japanese black pines",xlab="X",ylab="Y")




##############   Using the quadrat method to test for CSR
plot(simdat)
quadrat.test(simdat)
quadrat.test(simdat, 4, 3)

quadrat.test(simdat, alternative="regular")
quadrat.test(simdat, alternative="clustered")

# Using Monte Carlo p-values
quadrat.test(swedishpines) # Get warning, small expected values.
## Not run:
quadrat.test(swedishpines, method="M", nsim=4999)
quadrat.test(swedishpines, method="M", nsim=4999, conditional=FALSE)


################# K function

## Note that the Kest command computes an estimate of the K function and it allows for
## different correction methods: 
K.japanese <- Kest(japanesepines)
plot(K.japanese,main="K function \n Japanese black pines")

envjapK <- spatstat::envelope(japanesepines,fun=Kest,nrank=1,nsim=99)
plot(envjapK,main="K function \n Japanese black pines")


################# G function

## Note that the Gest command computes an estimate of the G function and it allows for
## different correction methods: 
G.japanese <- Gest(japanesepines)
plot(G.japanese,main="G function \n Japanese black pines") 
plot(spatstat::envelope(japanesepines,fun=Gest,nrank=1,nsim=99),col=rep(1,4),lwd=rep(2,4)) 



################# L function

L.japanese <- Lest(japanesepines)
plot(L.japanese,main="L function \n Japanese black pines") 
plot(spatstat::envelope(japanesepines,fun=Lest,nrank=1,nsim=99),col=rep(1,4),lwd=rep(2,4)) 



#################
#################


data(redwood)
summary(redwood)

## K function
K.redwood <- Kest(redwood)
plot(K.redwood,main="K function \n Redwood trees")

envredK <- envelope(redwood,fun=Kest,rank=1,nsim=99)
plot(envredK,main="K function \n Redwood trees")

## L function
plot(K.redwood, sqrt(iso/pi) ~ r, ylab="L(r)", main="L function \n Redwood trees")
abline(a=0,b=1,col="grey")


## Intensity function
int.redwood.01 <- density(redwood, 0.01)
plot(int.redwood.01,main="Kernel estimate of intensity function; h=0.01 \n Redwood seedlings and saplings",xlab="X",ylab="Y")

int.redwood.05 <- density(redwood, 0.05)
plot(int.redwood.05,main="Kernel estimate of intensity function; h=0.05 \n Redwood seedlings and saplings",xlab="X",ylab="Y")

int.redwood.1 <- density(redwood, 0.1)
plot(int.redwood.1,main="Kernel estimate of intensity function; h=0.1 \n Redwood seedlings and saplings",xlab="X",ylab="Y")

int.redwood.5 <- density(redwood, 0.5)
plot(int.redwood.5,main="Kernel estimate of intensity function; h=0.5 \n Redwood seedlings and saplings",xlab="X",ylab="Y")


#
redwood_point <- matrix(NA, nrow=length(redwood$x),ncol=2)
redwood_point[,1] <- redwood$x
redwood_point[,2] <- redwood$y

## Here we define a polygon that contains the events
redwood.poly <- as.points(c(min(redwood_point[,1]),max(redwood_point[,1]),max(redwood_point[,1]),min(redwood_point[,1])),
c(max(redwood_point[,2]),max(redwood_point[,2]),min(redwood_point[,2]),min(redwood_point[,2])))
plot(redwood.poly, type="l")


## This function computes the MSE for the estimated intensity function using different bandwidths
mse.redwood <- mse2d(redwood_point, poly=redwood.poly, nsmse=100, range=.5)

## Here we determine the optimal bandwidth as the one that minimizes the MSE
bw <- mse.redwood$h[which.min(mse.redwood$mse)] 
bw

plot(x=mse.redwood$h, y=mse.redwood$mse, xlab="Bandwidth", ylab="MSE", type="l",main="MSE as a function of bandwidth h \n Redwood seedlings and saplings")
i<-which.min(mse.redwood$mse)
points(mse.redwood$h[i], mse.redwood$mse[i])


## Here we plot the kernel estimate of the intensity function based on the selected bandwidth
int.redwood.bw <- density(redwood, bw)
plot(int.redwood.bw,main="Kernel estimate of intensity function; h=0.045 \n Redwood seedlings and saplings",xlab="X",ylab="Y")



###

data(cells)
summary(cells)

plot(cells$x,cells$y,type="p",pch=20,col="black",xlab="X",ylab="Y",main="Cell centers")
envcellsG <- envelope(cells,fun=Gest,rank=1,nsim=99)
plot(envcellsG,main="G function \n Cell centers")

envcellsF <- envelope(cells,fun=Fest,rank=1,nsim=99)
plot(envcellsF,main="F function \n Cell centers")

envcellsJ <- envelope(cells,fun=Jest,rank=1,nsim=99)
plot(envcellsJ,main="J function \n Cell centers")

K.cells <- Kest(cells)
plot(K.cells,main="K function \n Cell centers")

plot(K.cells, sqrt(iso/pi) ~ r, ylab="L(r)", main="L function \n Cell centers")
abline(a=0,b=1,col="grey")

## Intensity function
int.cells.01 <- density(cells, 0.01)
plot(int.cells.01,main="Kernel estimate of intensity function; h=0.01 \n Cell centers",xlab="X",ylab="Y")

int.cells.05 <- density(cells, 0.05)
plot(int.cells.05,main="Kernel estimate of intensity function; h=0.05 \n Cell centers",xlab="X",ylab="Y")

int.cells.1 <- density(cells, 0.1)
plot(int.cells.1,main="Kernel estimate of intensity function; h=0.1 \n Cell centers",xlab="X",ylab="Y")



cells_point <- matrix(NA, nrow=length(cells$x),ncol=2)
cells_point[,1] <- cells$x
cells_point[,2] <- cells$y

## Here we define a polygon that contains the events
cells.poly <- as.points(c(min(cells_point[,1]),max(cells_point[,1]),max(cells_point[,1]),min(cells_point[,1])),
c(max(cells_point[,2]),max(cells_point[,2]),min(cells_point[,2]),min(cells_point[,2])))
plot(cells.poly, type="l")


## This function computes the MSE for the estimated intensity function using different bandwidths
mse.cells <- mse2d(cells_point, poly=cells.poly, nsmse=100, range=.5)

## Here we determine the optimal bandwidth as the one that minimizes the MSE
bw <- mse.cells$h[which.min(mse.cells$mse)] 
bw

plot(x=mse.cells$h, y=mse.cells$mse, xlab="Bandwidth", ylab="MSE", type="l",main="MSE as a function of bandwidth h \n Cell centers")
i<-which.min(mse.cells$mse)
points(mse.cells$h[i], mse.cells$mse[i])


## Here we plot the kernel estimate of the intensity function based on the selected bandwidth
int.cells.bw <- density(cells, bw)
plot(int.cells.bw,main="Kernel estimate of intensity function; h=0.49 \n Cell centers",xlab="X",ylab="Y")


