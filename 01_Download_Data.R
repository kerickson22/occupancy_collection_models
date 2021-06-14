# Step 1: Download Raw Data

#Load packages

library(rgbif)
library(censusapi)

# Required files for this script:
# This script assumes you have downloaded the GADM
# shape file for the United States and stored it
# somewhere you can locate

#Download occurrence data

# 1. Tracheophytes WITH coordinates. #####
#This query
# was run on 13 February 2019 and includes records
# WITH coordinates that fall within the polygon
# `POLYGON((-88.90137 24.12911,-79.19824 24.12911,
# -79.19824 31.22069,-88.90137 31.22069,-88.90137
# 24.12911))` AND have no geospatial issues
#(`geospatial issues == FALSE`).

#The data is in the Darwin Core Archive Format
#The file size should be 64 MB
#https://doi.org/10.15468/dl.tzxjrm


res2 <- occ_download_get(key = "0040036-181108115102211", path = "./Raw Data/Tracheophytes with Coordinates", overwrite=T)
tracheophytes_coords <- occ_download_import(res2, fill=FALSE, quote="")
dim(tracheophytes_coords)
#There are only 214, 418 records
# (and the number of columns is correct)
save(tracheophytes_coords, file="./Raw Data/Tracheophytes with Coordinates/tracheophytes_coords.RData")


# 2. Tracheophytes WITHOUT coordinates #####
# This query was run on 12 February 2019 and
# includes records WITHOUT coordinates that have
# `stateProvince == florida`. The data is in
# the Darwin Core Archive Format
# The file size should be 116MB
# https://doi.org/10.15468/dl.zsgv2l

res2 <- occ_download_get(key = "0039523-181108115102211", path = "./Raw Data/Tracheophytes without Coordinates", overwrite=T)
tracheophytes_nocoords <- occ_download_import(res2, fill = FALSE, quote="")
dim(tracheophytes_nocoords)
#There are 294, 052 records

save(tracheophytes_nocoords, file="./Raw Data/Tracheophytes without Coordinates/tracheophytes_nocoords.RData")

#3. Shape file of the state of Florida with counties #####
# Load in a map of counties
# This step assumes you have downloaded the GADM
# shape file somewhere to your computer

counties <- rgdal::readOGR (
  './Raw Data/GADM/ver2pt8/WGS84',
  'USA_adm2'
)
florida <- counties[counties@data$NAME_1=="Florida", ]
save(florida, file="florida.RData")
rm(counties)

#Load in a state shapefile (for clipping points)
states <- rgdal::readOGR (
  './Raw Data/United States Shape',
  'USA_adm1'
)

florida_state <- states[states@data$NAME_1=="Florida", ]
save(florida_state, file="florida_state.RData")
rm(states)

#4. Census data for 2010 #####

#Prior to running this step, have to
# load census key into R environment
# (see details in the censusapi help manual
# for how to do this)

census_2010 <- getCensus(name="2010/dec/sf1?",
                         vars=c("NAME","COUNTY", "GEO_ID", "PCT001001"),
                         region="county:*", regionin = "state:12")
save(census_2010, file="census_2010.RData")








































