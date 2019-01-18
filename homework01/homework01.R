# Homework 01
setwd("~/WORKING_DIRECTORIES/biostat.696/homework01")
library(readr)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(akima)
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
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_gradient(name = "Observation Value") +
  coord_cartesian(xlim = c(-105,-65), ylim = c(25,50)) +
  theme_bw()

#
interp(long, lat, obs)