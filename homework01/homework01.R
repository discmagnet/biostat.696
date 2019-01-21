# Homework 01
setwd("~/WORKING_DIRECTORIES/biostat.696/homework01")
library(readr)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(akima)
library(fields)
data <- read_table2("data.txt")
colnames(data) <- c("long","lat","obs","CMAQ")

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