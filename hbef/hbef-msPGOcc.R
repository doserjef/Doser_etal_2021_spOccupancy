# hbef-msPGOcc.R: this script runs a nonspatial multispecies occupancy 
#                 model for the Hubbard Brook Experimental Forest foliage
#                 gleaning bird case study. 
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
#p.file.name <- "pfile-1"
p.file.name <- commandArgs(trailingOnly = TRUE)
chain <- unlist(strsplit(p.file.name, "-"))[2]
p.file <- read.table(paste("hbef/", p.file.name, sep = ''), 
		     sep = ' ', header = FALSE)
alpha.comm.start <- p.file[p.file[, 1] == 'alpha.comm', 2]
beta.comm.start <- p.file[p.file[, 1] == 'beta.comm', 2]
alpha.start <- p.file[p.file[, 1] == 'alpha', 2]
beta.start <- p.file[p.file[, 1] == 'beta', 2]
tau.beta.start <- p.file[p.file[, 1] == 'tau.beta', 2]
tau.alpha.start <- p.file[p.file[, 1] == 'tau.alpha', 2]

# Run Model ---------------------------------------------------------------
occ.ms.formula <- ~ scale(Elevation) + I(scale(Elevation)^2)
det.ms.formula <- ~ scale(day) + scale(tod) + I(scale(day)^2)
N <- dim(hbef2015$y)[1]
ms.inits <- list(alpha.comm = alpha.comm.start,
                    beta.comm = beta.comm.start,
                    beta = beta.start,
                    alpha = alpha.start,
                    tau.beta = tau.beta.start,
                    tau.alpha = tau.alpha.start,
                    z = apply(hbef2015$y, c(1, 2), max, na.rm = TRUE))
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72), 
                  alpha.comm.normal = list(mean = 0, var = 2.72),
                  tau.beta.ig = list(a = 0.1, b = 0.1),
                  tau.alpha.ig = list(a = 0.1, b = 0.1)) 
n.samples <- 150000
n.burn <- 50000
n.thin <- 40
# Run the model
out <- msPGOcc(occ.formula = occ.ms.formula, 
	       det.formula = det.ms.formula, 
	       data = hbef2015, 
	       inits = ms.inits, 
	       n.samples = n.samples,
	       priors = ms.priors, 
	       n.omp.threads = 1, 
	       verbose = TRUE, 
	       n.burn = n.burn, 
	       n.thin = n.thin, 
	       n.report = 1500)
n.post <- length(seq(from = n.burn + 1, 
		     to = n.samples, 
		     by = as.integer(n.thin)))
# Save the results
save(out, file = paste("results/hbef-msPGOcc-", chain, "-", n.post, "-", 
		       Sys.Date(), ".R", sep = ''))

