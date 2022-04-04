# bbs-PGOcc.R: this script runs a single species occupancy model for 
#              Black-throated Blue Warbler across the eastern US in 2018. 
# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(spOccupancy)
library(coda)

# Read in the data --------------------------------------------------------
load("bbs/data/bbs-btnw-bundle.R")

# Get inits values -----------------------------------------------------
# This is used to run scripts from the command line to run multiple chains
# simultaneously across multiple cores. This is a basic way of running 
# multiple chains in parallel using spOccupancy. To run a single chain 
# directly in the script, uncomment the following line of code and comment
# out the line of code underneath it: 
#p.file.name <- "pfile-1"
p.file.name <- commandArgs(trailingOnly = TRUE)
chain <- unlist(strsplit(p.file.name, "-"))[2]
p.file <- read.table(paste("bbs/", p.file.name, sep = ''),
		     sep = ' ', header = FALSE)
alpha.start <- p.file[p.file[, 1] == 'alpha', 2]
beta.start <- p.file[p.file[, 1] == 'beta', 2]

# Run Model ---------------------------------------------------------------
occ.formula <- ~ elev + elev.2 + pf
det.formula <- ~ day + day.2 + tod + (1 | obs)
p.det <- length(bbs.btnw.dat$det.covs)
p.occ <- ncol(bbs.btnw.dat$occ.covs) + 1
inits <- list(alpha = rep(alpha.start, p.det),
                 beta = rep(beta.start, p.occ),
                 z = apply(bbs.btnw.dat$y, 1, max, na.rm = TRUE))
priors <- list(beta.normal = list(mean = rep(0, p.occ),
                                  var = rep(2.72, p.occ)),
               alpha.normal = list(mean = rep(0, p.det),
                                   var = rep(2.72, p.det)))
n.samples <- 50000
n.burn <- 10000
n.thin <- 20

# Run the model
out <- PGOcc(occ.formula = occ.formula,
	     det.formula = det.formula,
	     data = bbs.btnw.dat,
	     inits = inits,
	     n.samples = n.samples,
	     priors = priors,
	     n.omp.threads = 1,
	     verbose = TRUE,
	     n.burn = n.burn,
	     n.thin = n.thin,
	     n.report = 500)
# Save results
save(out, file = paste("results/bbs-PGOcc-", chain, "-", 
		       Sys.Date(), ".R", sep = ''))

