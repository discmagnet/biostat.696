# Homework 01
setwd("~/WORKING_DIRECTORIES/biostat.696/homework01")
library(readr)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(akima)
library(fields)
library(dplyr)
library(broom)
library(stargazer)
data <- read_table2("data.txt")
colnames(data) <- c("long","lat","obs","CMAQ")

R <- 6371
data <- mutate(data, X = R*cos(lat)*cos(long),
                     Y = R*cos(lat)*sin(long),
                     Z = R*sin(lat))

# Display the monitoring locations and corresponding
# observed values in the United States.
usa <- map_data("usa")
states <- map_data("state")
ggplot() +
  geom_polygon(data = states,
               aes(x=long, y=lat, group=group),
               fill = "white", color = "gray") +
  geom_polygon(data = usa,
               aes(x=long, y=lat, group=group),
               fill = NA, color = "black") +
  geom_point(data = data,
             aes(x=long, y = lat, color = obs)) +
  ggtitle("Ozone Concentrations in the United States") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_gradient(name = "Observation Value") +
  coord_cartesian(xlim = c(-105,-65), ylim = c(25,50)) 
  

# Plot the interpolated surface for the observed ozone
# concentration and add contour lines to the plot.
source("heatmap.R")
surf.ozone <- heatmap(data$long, data$lat, data$obs)
ggplot(surf.ozone, aes(x=lon,y=lat)) +
  geom_polygon(data = states,
               aes(x=long, y=lat, group=group),
               fill = "white", color = "black") +
  geom_raster(aes(fill = obs), interpolate = T) +
  geom_contour(aes(z = obs), color = "black") +
  geom_polygon(data = usa,
               aes(x=long, y=lat, group=group),
               fill = NA, color = "black") +
  geom_polygon(data = states,
               aes(x=long, y=lat, group=group),
               fill = NA, color = "black") +
  ggtitle("Ozone Concentrations in the United States") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_cartesian(xlim = c(-105,-65), ylim = c(25,50)) +
  scale_fill_gradient(na.value = NA, low = "yellow", high = "red", name = "Observation Value")

# EDA of Ozone Data
# Histogram of observed ozone concentration
hist(data$obs, xlab = "Ozone Concentration", ylab = "Frequency",
     main = "Histogram of Observed Ozone Concentration")

data$obs_sqrt <- sqrt(data$obs)
hist(data$obs_sqrt, xlab = "Square Root of Ozone Concentration",
     ylab = "Frequency", main = "Histogram of Transformed Observed Ozone Concentration")

data$obs_log <- log(data$obs)
hist(data$obs_log, xlab = "Log of Ozone Concentration",
     ylab = "Frequency", main = "Histogram of Transformed Observed Ozone Concentration")

# Linear Regression of Ozone Data
# consider a cubic trend surface model
data$X2 <- (data$X)^2
data$Y2 <- (data$Y)^2
data$XY <- (data$X)*(data$Y)
data$XZ <- (data$X)*(data$Z)
data$YZ <- (data$Y)*(data$Z)
data$X2Y <- ((data$X)^2)*(data$Y)
data$XY2 <- (data$X)*((data$Y)^2)
data$X2Z <- ((data$X)^2)*(data$Z)
data$Y2Z <- ((data$Y)^2)*(data$Z)
data$X3 <- (data$X)^3
data$Y3 <- (data$Y)^3

mod01 <- lm(data = data,
            obs ~ CMAQ+X+Y+Z+X2+Y2+XY+XZ+YZ+X2Y+XY2+X2Z+Y2Z+X3+Y3)
summary(mod01)
mod01_coeff <- summary(mod01)$coeff[1:14]
data$res <- mod01$residuals

beta <- tidy(mod01)
beta$p.value <- round(beta$p.value,3)
beta$p.value[beta$p.value==0] <- "<0.001"

# Semi-variogram
# empirical variogram
emp_var <- variogram(res~1, locations =~ X+Y+Z, data, width = 50, cutoff = 6000)
emp_var
plot(emp_var,col="black",type="p",pch=20, main = "Empirical Variogram")
# empirical directional variogram
dir_var <- variogram(res~1,locations =~ X+Y+Z, data, width = 50, cutoff = 6000,
                     alpha=c(0,45,90,135))
dir.E <- dir_var[dir_var$dir.hor==0,]
dir.NE <- dir_var[dir_var$dir.hor==45,]
dir.N <- dir_var[dir_var$dir.hor==90,]
dir.NW <- dir_var[dir_var$dir.hor==135,]
plot(dir.E$dist,dir.E$gamma,col="black",type="p",pch=20,xlab="Distance",ylab="Semivariance",
     main = "E-W Empirical Directional Variogram", ylim=c(0,300))
plot(dir.NE$dist,dir.NE$gamma,col="black",type="p",pch=20,xlab="Distance",ylab="Semivariance",
     main = "NE-SW Empirical Directional Variogram", ylim=c(0,300))
plot(dir.N$dist,dir.N$gamma,col="black",type="p",pch=20,xlab="Distance",ylab="Semivariance",
     main = "N-S Empirical Directional Variogram", ylim=c(0,300))
plot(dir.NW$dist,dir.NW$gamma,col="black",type="p",pch=20,xlab="Distance",ylab="Semivariance",
     main = "NW-SE Empirical Directional Variogram", ylim=c(0,300))

exp_var <- fit.variogram(emp_var,vgm(psill=50,"Mat",1000,100),fit.method=2)
exp_var
plot(emp_var, exp_var)
