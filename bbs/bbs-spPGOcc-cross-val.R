# bbs-spPGOcc-cross-val.R: this script runs a single species spatial occupancy 
#                          model with 10-fold cross-validation for Black-throated
#                          Blue Warbler across the eastern US in 2018. The model 
#                          is fit with a nearest neighbor Gaussian process.
# Author: Jeffrey W. Doser
# Citation: 
rm(list = ls())
library(spOccupancy)
library(coda)
library(sf)

# Read in the data --------------------------------------------------------
load("bbs/data/bbs-btnw-bundle.R")
 
# Get coordinates in a projection rather than lat-long. 
coords.sf <- st_as_sf(data.frame(bbs.btnw.dat$coords), 
		      coords = c("Longitude", "Latitude"), 
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Albers equal area across contiguous US. 
coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Get coordinates in Albers Equal Area
coords.albers <- st_coordinates(coords.sf.albers)
# Convert coordinates to km in Albers equal area.
bbs.btnw.dat$coords <- coords.albers / 1000

# Get inits values -----------------------------------------------------
p.file.name <- "pfile-sp-3"
p.file <- read.table(paste("bbs/", p.file.name, sep = ''),
		     sep = ' ', header = FALSE)
alpha.start <- p.file[p.file[, 1] == 'alpha', 2]
beta.start <- p.file[p.file[, 1] == 'beta', 2]
sigma.sq.start <- p.file[p.file[, 1] == 'sigma.sq', 2]
phi.start <- p.file[p.file[, 1] == 'phi', 2]
w.start <- p.file[p.file[, 1] == 'w', 2]

# Run Model ---------------------------------------------------------------
occ.formula <- ~ elev + elev.2 + pf
det.formula <- ~ day + day.2 + tod + (1 | obs)
p.det <- length(bbs.btnw.dat$det.covs)
p.occ <- ncol(bbs.btnw.dat$occ.covs) + 1
# Prep spatial stuff
dist.bbs <- dist(bbs.btnw.dat$coords)
mean.dist <- mean(dist.bbs)
min.dist <- min(dist.bbs)
max.dist <- max(dist.bbs)
inits <- list(alpha = rep(alpha.start, p.det),
                 beta = rep(beta.start, p.occ),
		 sigma.sq = sigma.sq.start, 
		 phi = phi.start, 
		 w = rep(w.start, nrow(bbs.btnw.dat$y)),
                 z = apply(bbs.btnw.dat$y, 1, max, na.rm = TRUE))
priors <- list(beta.normal = list(mean = rep(0, p.occ),
                                  var = rep(2.72, p.occ)),
               alpha.normal = list(mean = rep(0, p.det),
                                   var = rep(2.72, p.det)), 
	       phi.unif = c(3 / max.dist, 3 / min.dist), 
	       sigma.sq.ig = c(2, 5))
batch.length <- 25
n.batch <- 2000
n.burn <- 10000
n.thin <- 20
n.report <- 20
tuning <- list(phi = 1)
# Run the model
out <- spPGOcc(occ.formula = occ.formula,
	       det.formula = det.formula,
	       data = bbs.btnw.dat,
	       inits = inits,
	       batch.length = batch.length,
	       n.batch = n.batch, 
	       tuning = tuning,
	       priors = priors,
	       n.omp.threads = 5,
	       verbose = TRUE,
	       NNGP = TRUE, 
	       cov.model = 'exponential', 
	       n.neighbors = 10,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.report = n.report, 
	       k.fold = 10, 
	       k.fold.threads = 10)
# Save the results
save(out, file = paste("results/bbs-spPGOcc-cross-val-", 
		       Sys.Date(), ".R", sep = ''))

