# Homework 04 R Code
setwd("~/WORKING_DIRECTORIES/biostat.696/homework04")

# Load Libraries
library(readr)
library(spatstat)
library(dplyr)
library(smacpod)

# Import Dataset
cancer <- read_table2("Cancer_data.txt")
colnames(cancer) <- c("x","y","type") # Rename columns
cancer$type <- factor(cancer$type) # Factor w/ 2 levels: "larynx" "lung"

# Import Polygon to Define Window of Interest
polygon <- read_table2("Polygon_cancer_data.txt")
colnames(polygon) <- c("x","y")

# Locations of the larynx cases
larynx <- subset(cancer, type == levels(type)[1])
larynx$xy <- paste(larynx$x,larynx$y)
larynx <- distinct(larynx,larynx$xy,.keep_all = T)[,1:2]

# Locations of the lung cases
lung <- subset(cancer, type == levels(type)[2])
lung$xy <- paste(lung$x,lung$y)
lung <- distinct(lung,lung$xy,.keep_all = T)[,1:2]

# Plot the Larynx Cases
plot(larynx, main = "Larynx Cancer Cases in the UK",
     xlab = "Easting (km)",
     ylab = "Northing (km)",
     asp = 1)
polygon(polygon, col = "gray")
points(larynx[,1:2], pch = 20)

# Set up ppp objects
poly <- owin(poly = polygon)
larynx.ppp <- as.ppp(larynx, W = poly)
lung.ppp <- as.ppp(lung, W = poly)
# ppp object used in part (h)
all <- rbind(larynx, lung)
all.ppp <- as.ppp(all, W = poly)
all.ppp$markformat <- "vector"
all.ppp$marks <- factor(c(rep("larynx",nrow(larynx)),rep("lung",nrow(lung))))

# Estimate the K function
K.larynx <- Kest(larynx.ppp)
plot(K.larynx, main = "K function \n Larynx Cancer Cases")
env_larynx_K <- envelope(larynx.ppp, fun = Kest, nrank = 1, nsim = 99)
plot(env_larynx_K, main = "K function \n Larynx Cancer Cases")

# Estimate the L function
L.larynx <- Lest(larynx.ppp)
plot(L.larynx,main="L function \n Larynx Cancer Cases")
env_larynx_L <- envelope(larynx.ppp, fun = Lest, nrank=1, nsim=99)
plot(env_larynx_L, main = "L function \n Larynx Cancer Cases", col=rep(1,4), lwd=rep(2,4))

# Test whether the point pattern is a Poisson point process
quadrat.test(larynx.ppp)
quadrat.test(larynx.ppp, alternative = "regular")
quadrat.test(larynx.ppp, alternative = "clustered")
quadrat.test(larynx.ppp, method="MonteCarlo", nsim=4999)

# Estimate the intensity function for larynx cases
int.larynx <- density(larynx.ppp, bw = "bw.diggle")
plot(int.larynx, main = "Kernel estimate of intensity function \n Larynx Cancer Cases")

# Plot the Lung Cases
plot(lung, main = "Lung Cancer Cases in the UK",
     xlab = "Easting (km)",
     ylab = "Northing (km)",
     asp = 1)
polygon(polygon, col = "gray")
points(larynx[,1:2], pch = 20)

# Estimate the intensity function for lung cases (controls)
int.lung <- density(lung.ppp, bw = "bw.diggle")
plot(int.lung, main = "Kernel estimate of intensity function \n Lung Cancer Cases")

# Fit a log-Gaussian Cox process
larynx.logcp <- lgcp.estK(larynx.ppp, c(sigma2=10, alpha=2))
larynx.logcp
plot(larynx.logcp,main="Fitted K function and theoretical K function \n Larynx Cancer Cases, log Gaussian-Cox process")

# Fit a Matern cluster point process
larynx.matclust <- matclust.estK(larynx.ppp, c(kappa=0.1, R=2))
larynx.matclust
plot(larynx.matclust,main="Fitted K function and theoretical K function \n Larynx Cancer Cases, Matern cluster process")

# Perform a Kutzick-Edwards test
qnn.test(all.ppp, q = c(3,5,7,9), case = 1)