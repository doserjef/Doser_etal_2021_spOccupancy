# hbef-msPGOcc-int.R: this script runs a nonspatial multispecies occupancy 
#                     model with an intercept only occurrence model for the 
#                     Hubbard Brook Experimental Forest foliage-gleaning 
#                     bird case study. 
# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(spOccupancy)
library(coda)

# Load data ---------------------------------------------------------------
data(hbef2015)

# Set inits values -----------------------------------------------------
alpha.comm.start <- 0
beta.comm.start <- 0
alpha.start <- 0
beta.start <- 0
tau.sq.beta.start <- 1
tau.sq.alpha.start <- 1

# Prep for Spatial Models -------------------------------------------------
# Distances between sites
dist.hbef <- dist(hbef2015$coords)
min.dist <- min(dist.hbef)
max.dist <- max(dist.hbef)

# Run Model ---------------------------------------------------------------
# Intercept Only ----------------------
N <- dim(hbef2015$y)[1]
ms.inits <- list(alpha.comm = alpha.comm.start,
                    beta.comm = beta.comm.start,
                    beta = beta.start,
                    alpha = alpha.start,
                    tau.sq.sq.beta = tau.sq.beta.start,
                    tau.sq.sq.alpha = tau.sq.alpha.start,
                    z = apply(hbef2015$y, c(1, 2), max, na.rm = TRUE))
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72), 
                  alpha.comm.normal = list(mean = 0, var = 2.72),
                  tau.sq.sq.beta.ig = list(a = 0.1, b = 0.1),
                  tau.sq.sq.alpha.ig = list(a = 0.1, b = 0.1)) 

n.samples <- 150000
n.burn <- 50000
n.thin <- 100
out.int <- msPGOcc(occ.formula = ~ 1, 
	           det.formula = ~ scale(day) + scale(tod) + I(scale(day)^2), 
	           data = hbef2015, 
	           inits = ms.inits, 
	           n.samples = n.samples,
	           priors = ms.priors, 
	           n.omp.threads = 1, 
	           verbose = TRUE, 
	           n.burn = n.burn, 
	           n.thin = n.thin, 
	           n.report = 1500, 
                   k.fold = 4, 
                   k.fold.threads = 4) 

save(out.int, file = paste("results/hbef-msPGOcc-int", "-", 
		            Sys.Date(), ".R", sep = ''))

