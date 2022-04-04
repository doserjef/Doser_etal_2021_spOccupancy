# bbs-data-prep.R: this script takes the raw BBS data and converts it into
#                  the necessary form for analysis in spOccupancy for an 
#                  analysis of Black-throated Green Warbler across the 
#                  eastern US. 
# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(tidyverse)
library(lubridate)
library(sp)
library(raster)
library(FedData)
library(sf)
library(stars)

# Read in BBS Data --------------------------------------------------------
# These are the 10-stop summary data for all states included in analysis.
bbs.dat <- list.files(path = "bbs/data/btbw-states/", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()
# Get associated route data
route.dat <- read.csv("bbs/data/routes.csv")
# Get associated weather data
weather.dat <- read.csv("bbs/data/weather.csv")
# Note that route id is nested within state. 
# Join BBS data with route data
bbs.dat <- left_join(bbs.dat, route.dat, by = c('Route', 'CountryNum', 'StateNum'))
# Only grab data from 2018
bbs.2018 <- bbs.dat %>%
  filter(Year == 2018)
# Select columns of interest
bbs.2018 <- bbs.2018 %>%
  dplyr::select(RouteDataID, Latitude, Longitude, AOU, starts_with("Count")) %>%
  dplyr::select(-CountryNum)
# Fill in implicit zeros
bbs.2018 <- bbs.2018 %>%
  complete(AOU, nesting(RouteDataID, Latitude, Longitude))
# Replace NAs with 0s for all columns at once.
bbs.2018 <- bbs.2018 %>%
  replace(is.na(.), 0)
# Only grab species of interest
# Code for BTNW is 6670
curr.sp <- 6670
bbs.2018 <- bbs.2018 %>%
  filter(AOU == curr.sp)

# Format Detection and raw data -------------------------------------------
# All detection level variables are at the stop level.
# Join data with Weather data
# Get date in proper format
weather.dat <- weather.dat %>% 
  unite('date', sep = '-', Year, Month, Day, remove = FALSE)
weather.dat$date <- as.Date(weather.dat$date, tz = "America/New_York")
# Get julian date of each survey
weather.dat$julian <- as.numeric(format(weather.dat$date, '%j'))
btnw.dat <- left_join(bbs.2018, weather.dat, by = c('RouteDataID'))

# Create arrays for packaging in a list bundle for spOccupancy. 
# Detection-nondetection data
y <- as.matrix(btnw.dat[, c('Count10', 'Count20', 'Count30', 'Count40', 'Count50')])
y <- ifelse(y > 0, 1, 0)
# Detection covariates
det.covs <- list(day = c(scale(btnw.dat$julian)), 
		 day.2 = c(scale(btnw.dat$julian)^2), 
		 tod = c(scale(btnw.dat$StartTime)), 
		 obs = as.numeric(factor(btnw.dat$ObsN)))


# Occurrence covariates ---------------------------------------------------
coords <- btnw.dat[, c('Latitude', 'Longitude')]
# The code below takes a while to run, and is probably best done on a remote
# machine if possible. Need to break this up into pieces to make it work, 
# otherwise it times out
coords.sp <- data.frame(coords, val = apply(y, 1, max))
coords.sp <- coords.sp %>% arrange(Longitude, Latitude)
coordinates(coords.sp) <- ~Longitude + Latitude
proj4string(coords.sp) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
# Loop through all the sites
J <- nrow(coords)
vals <- split(1:J, ceiling(seq_along(1:J)/4))
elev.btnw <- rep(0, J)
for.btnw <- rep(0, J)
# Get proportion of forest
props <- function(a, na.rm = TRUE) {
  my.sum <- sum(!is.na(a))	
  prop.for <- sum(a %in% c(41, 42, 43), na.rm = na.rm) / my.sum
  return(prop.for)
}
ned.dat <- list()
nlcd.dat <- list()
for (i in 1:length(vals)) {
  print(paste("Currently on iteration ", i, " out of ", length(vals), sep = ''))
  ned.dat[[i]] <- get_ned(template = coords.sp[vals[[i]], ], label = paste('btnw', i))
  nlcd.dat[[i]] <- get_nlcd(template = coords.sp[vals[[i]], ], label = paste('btnw', i), year = 2016)
  elev.btnw[vals[[i]]] <- extract(ned.dat[[i]], coords.sp[vals[[i]], ])
  for.btnw[vals[[i]]] <- extract(nlcd.dat[[i]], coords.sp[vals[[i]], ], buffer = 1000, fun = props)
}

# Get everything in the same order -----------------------------------------
all.vars <- data.frame(day = det.covs$day, 
		       day.2 = det.covs$day.2, 
		       tod = det.covs$tod, 
		       obs = det.covs$obs, 
		       y, 
		       coords) 
# Temporary
tmp.vars <- data.frame(y, coords)
tmp.vars <- tmp.vars %>%
  arrange(Longitude, Latitude)
y <- as.matrix(tmp.vars[, c('Count10', 'Count20', 'Count30', 'Count40', 'Count50')])

all.vars <- all.vars %>%
  arrange(Longitude, Latitude)

det.covs <- list(day = all.vars$day, 
		 day.2 = all.vars$day.2, 
		 tod = all.vars$tod, 
		 obs = all.vars$obs)
y <- as.matrix(all.vars[, c('Count10', 'Count20', 'Count30', 'Count40', 'Count50')])

coords <- as.matrix(all.vars[, c("Longitude", "Latitude")])
elev.mean <- mean(elev.btnw, na.rm = TRUE)
elev.sd <- sd(elev.btnw, na.rm = TRUE)
for.mean <- mean(for.btnw, na.rm = TRUE)
for.sd <- sd(for.btnw, na.rm = TRUE)

occ.covs <- data.frame(elev = c(scale(elev.btnw)), 
		       pf = c(scale(for.btnw)), 
		       elev.2 = c(scale(elev.btnw))^2)

bbs.btnw.dat <- list(y = y, det.covs = det.covs, 
		     occ.covs = occ.covs, coords = coords)
# Save data in a list bundle for spOccupancy.
save(bbs.btnw.dat, file = "data/bbs-btnw-bundle.R")
