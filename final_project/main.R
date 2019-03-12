# Final Project
setwd("~/WORKING_DIRECTORIES/biostat.696/final_project")

# Import Data and Export a single RData file
library(readr)
storm96 <- read_csv("stormdata/Stormdata_1996.csv")
storm97 <- read_csv("stormdata/Stormdata_1997.csv")
storm98 <- read_csv("stormdata/Stormdata_1998.csv")
storm99 <- read_csv("stormdata/Stormdata_1999.csv")
storm00 <- read_csv("stormdata/Stormdata_2000.csv")
storm01 <- read_csv("stormdata/Stormdata_2001.csv")
storm02 <- read_csv("stormdata/Stormdata_2002.csv")
storm03 <- read_csv("stormdata/Stormdata_2003.csv")
storm04 <- read_csv("stormdata/Stormdata_2004.csv")
storm05 <- read_csv("stormdata/Stormdata_2005.csv")
storm06 <- read_csv("stormdata/Stormdata_2006.csv")
storm07 <- read_csv("stormdata/Stormdata_2007.csv")
storm08 <- read_csv("stormdata/Stormdata_2008.csv")
storm09 <- read_csv("stormdata/Stormdata_2009.csv")
storm10 <- read_csv("stormdata/Stormdata_2010.csv")
storm11 <- read_csv("stormdata/Stormdata_2011.csv")
storm12 <- read_csv("stormdata/Stormdata_2012.csv")
storm13 <- read_csv("stormdata/Stormdata_2013.csv")
storm_all <- rbind(storm96,storm97,storm98,storm99,storm00,
               storm01,storm02,storm03,storm04,storm05,
               storm06,storm07,storm08,storm09,storm10,
               storm11,storm12,storm13)
tornado <- subset(storm_all, EVENT_TYPE == "Tornado")
final_data <- tornado[c("BEGIN_DATE_TIME","END_DATE_TIME",
                        "EVENT_ID","STATE","INJURIES_DIRECT",
                        "DEATHS_DIRECT","DAMAGE_PROPERTY",
                        "DAMAGE_CROPS","TOR_F_SCALE","TOR_LENGTH",
                        "TOR_WIDTH","BEGIN_LAT","BEGIN_LON",
                        "END_LAT","END_LON")]
# Remove weak tornados (Fujita Scale is F0 or EF0 or NA) - 18,458
final_data2 <- final_data[!is.na(final_data$TOR_F_SCALE),]
final_data3 <- final_data2[final_data2$TOR_F_SCALE != "F0",]
final_data4 <- final_data3[final_data3$TOR_F_SCALE != "EF0",]
# Remove tornados without coordinates - 124
final_data5 <- final_data4[!is.na(final_data4$END_LAT),]
# Note that there are still missing values in this dataset,
# just not sure which cases or variables we want to exclude
# 27 - missing TOR_LENGTH and TOR_WIDTH
# 1,072 - missing DAMAGE_PROPERTY
# 3,374 - missing DAMAGE_CROPS

tordata <- final_data5

# Clean-up variables
rm(storm00,storm01,storm02,storm03,storm04,storm05,storm06,storm07,
   storm08,storm09,storm10,storm11,storm12,storm13,storm99,storm98,
   storm97,storm96,storm_all,tornado,final_data,final_data2,
   final_data3,final_data4,final_data5)

save.image("~/WORKING_DIRECTORIES/biostat.696/final_project/tordata.RData")