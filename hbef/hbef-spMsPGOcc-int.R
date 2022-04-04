# hbef-spMsPGOcc-int.R: this script runs a spatial multispecies occupancy 
#                       model with an intercept only occurrence model for the 
#                       Hubbard Brook Experimental Forest foliage-gleaning 
#                       bird case study. 
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
# Spatial Intercept Only ---------------
N <- dim(hbef2015$y)[1]
phi.start <- 3 / mean(dist.hbef)
sigma.sq.start <- 1
w.start <- 0
cov.model <- "exponential"
ms.inits <- list(alpha.comm = alpha.comm.start, 
                    beta.comm = beta.comm.start, 
                    beta = beta.start, 
                    alpha = alpha.start,
                    tau.sq.sq.beta = tau.sq.beta.start, 
                    tau.sq.sq.alpha = tau.sq.alpha.start, 
                    z = apply(hbef2015$y, c(1, 2), max, na.rm = TRUE), 
		    sigma.sq = sigma.sq.start, 
		    phi = phi.start, 
		    w = matrix(w.start, N, dim(hbef2015$y)[2]))
min.dist <- min(dist.hbef)
max.dist <- max(dist.hbef)
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72),
                  tau.sq.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.sq.alpha.ig = list(a = 0.1, b = 0.1), 
		  sigma.sq.ig = list(a = 2, b = 1),
		  phi.unif = list(a = 3/max.dist, b = 3 / min.dist))
batch.length <- 25
n.batch <- 6000
n.burn <- 50000
n.thin <- 100
ms.tuning <- list(phi = 1)
out.sp.int <- spMsPGOcc(occ.formula = ~ 1, 
		        det.formula = ~ scale(day) + scale(tod) + I(scale(day)^2),
		        data = hbef2015, 
		        inits = ms.inits, 
		        n.batch = n.batch, 
		        batch.length = batch.length, 
		        accept.rate = 0.43, 
		        priors = ms.priors, 
		        cov.model = cov.model, 
		        tuning = ms.tuning, 
		        n.omp.threads = 5, 
		        verbose = TRUE, 
		        NNGP = TRUE, 
		        n.neighbors = 5,
		        n.report = 60, 
		        n.burn = n.burn, 
		        n.thin = n.thin) 

save(out.sp.int, file = paste("results/hbef-spMsPGOcc-int-no-cross", "-", 
		               Sys.Date(), ".R", sep = ''))
