# hbef-spMPGOcc.R: this script runs a spatial multispecies occupancy 
#                  model for the Hubbard Brook Experimental Forest foliage
#                  gleaning bird case study. 
# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(spOccupancy)
library(coda)

# Load data ---------------------------------------------------------------
data(hbef2015)

# Get inits values -----------------------------------------------------
# This is used to run scripts from the command line to run multiple chains
# simultaneously across multiple cores. This is a basic way of running 
# multiple chains in parallel using spOccupancy. To run a single chain 
# directly in the script, uncomment the following line of code and comment
# out the line of code underneath it: 
p.file.name <- "pfile-sp-1"
#p.file.name <- commandArgs(trailingOnly = TRUE)
chain <- unlist(strsplit(p.file.name, "-"))[3]
# This should be user input
p.file <- read.table(paste("hbef/", p.file.name, sep = ''), 
		     sep = ' ', header = FALSE)
alpha.comm.start <- p.file[p.file[, 1] == 'alpha.comm', 2]
beta.comm.start <- p.file[p.file[, 1] == 'beta.comm', 2]
alpha.start <- p.file[p.file[, 1] == 'alpha', 2]
beta.start <- p.file[p.file[, 1] == 'beta', 2]
sigma.sq.start <- p.file[p.file[, 1] == 'sigma.sq', 2]
phi.start <- p.file[p.file[, 1] == 'phi', 2]
tau.beta.start <- p.file[p.file[, 1] == 'tau.beta', 2]
tau.alpha.start <- p.file[p.file[, 1] == 'tau.alpha', 2]
w.start <- p.file[p.file[, 1] == 'w', 2]

# Run Model ---------------------------------------------------------------
occ.formula <- ~ scale(Elevation) + I(scale(Elevation)^2)
det.formula <- ~ scale(day) + scale(tod) + I(scale(day)^2)
# Number of species
N <- dim(hbef2015$y)[1]
# Distances between sites
dist.hbef <- dist(hbef2015$coords)
# Exponential covariance model
cov.model <- "exponential"
ms.inits <- list(alpha.comm = alpha.comm.start, 
                    beta.comm = beta.comm.start, 
                    beta = beta.start, 
                    alpha = alpha.start,
                    tau.beta = tau.beta.start, 
                    tau.alpha = tau.alpha.start, 
                    z = apply(hbef2015$y, c(1, 2), max, na.rm = TRUE), 
		    sigma.sq = sigma.sq.start, 
		    phi = phi.start, 
		    w = matrix(w.start, N, dim(hbef2015$y)[2]))
min.dist <- min(dist.hbef)
max.dist <- max(dist.hbef)
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72),
                  tau.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.alpha.ig = list(a = 0.1, b = 0.1), 
		  sigma.sq.ig = list(a = 2, b = 1),
		  phi.unif = list(a = 3/max.dist, b = 3 / min.dist))
batch.length <- 25
n.batch <- 6000
n.burn <- 50000
n.thin <- 40
n.report <- 60
ms.tuning <- list(phi = 1)
# Values for reporting
# Run model
out <- spMsPGOcc(occ.formula = occ.formula, 
		 det.formula = det.formula, 
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
		 n.report = n.report, 
		 n.burn = n.burn, 
		 n.thin = n.thin)
# Read results
save(out, file = paste("results/hbef-spMsPGOcc-", chain, "-", 
		       Sys.Date(), ".R", sep = ''))

