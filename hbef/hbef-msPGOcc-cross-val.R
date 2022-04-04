# hbef-msPGOcc-cross-val.R: this script runs a nonspatial multispecies occupancy 
#                           model for the Hubbard Brook Experimental Forest foliage
#                           gleaning bird case study. The model is run using a four-fold
#                           cross-validation procedure for model comparison. 
# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(spOccupancy)
library(coda)

# Load data ---------------------------------------------------------------
data(hbef2015)

# Get inits values -----------------------------------------------------
p.file.name <- "pfile-3"
chain <- unlist(strsplit(p.file.name, "-"))[2]
# This should be user input
p.file <- read.table(paste("hbef/", p.file.name, sep = ''), 
		     sep = ' ', header = FALSE)
alpha.comm.start <- p.file[p.file[, 1] == 'alpha.comm', 2]
beta.comm.start <- p.file[p.file[, 1] == 'beta.comm', 2]
alpha.start <- p.file[p.file[, 1] == 'alpha', 2]
beta.start <- p.file[p.file[, 1] == 'beta', 2]
tau.sq.beta.start <- p.file[p.file[, 1] == 'tau.sq.beta', 2]
tau.sq.alpha.start <- p.file[p.file[, 1] == 'tau.sq.alpha', 2]

# Run Model ---------------------------------------------------------------
occ.ms.formula <- ~ scale(Elevation) + I(scale(Elevation)^2)
det.ms.formula <- ~ scale(day) + scale(tod) + I(scale(day)^2)
N <- dim(hbef2015$y)[1]
ms.inits <- list(alpha.comm = alpha.comm.start,
                    beta.comm = beta.comm.start,
                    beta = beta.start,
                    alpha = alpha.start,
                    tau.sq.beta = tau.sq.beta.start,
                    tau.sq.alpha = tau.sq.alpha.start,
                    z = apply(hbef2015$y, c(1, 2), max, na.rm = TRUE))
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72), 
                  alpha.comm.normal = list(mean = 0, var = 2.72),
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1)) 
n.samples <- 150000
n.burn <- 50000
n.thin <- 100
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
	       n.report = 1500, 
	       k.fold = 4, 
	       k.fold.threads = 4)
n.post <- length(seq(from = n.burn + 1, 
		     to = n.samples, 
		     by = as.integer(n.thin)))


save(out, file = paste("results/hbef-msPGOcc-cross-val-", chain, "-", n.post, "-", 
		       Sys.Date(), ".R", sep = ''))

