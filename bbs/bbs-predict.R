# bbs-predict.R: this script predicts occurrence for the Black-throated
#                Green Warbler across the eastern US using a spatial 
#                occupancy model with an NNGP. 
rm(list = ls())
library(spOccupancy)
library(coda)

# Read in spatial model
load("results/bbs-spPGOcc-3-2022-03-22.R")

# Read in prediction values
load("bbs/data/bbs-pred-data.rda")
coords.0 <- pred.df[, c('X', 'Y')]
# Standardize by values used to fit the model
elev.pred <- (pred.df$elev - 261.794) / 182.3823
for.pred <- (pred.df$forest - 0.4021555) / 0.2991306
X.0 <- cbind(1, elev.pred, elev.pred^2, for.pred)
names(X.0) <- c('int', 'elev', 'elev.2', 'pf')
# Change the number of threads accordingly.
out.pred <- predict(out, X.0, coords.0, n.omp.threads = 10)

psi.0.samples <- out.pred$psi.0.samples
save(psi.0.samples, coords.0, file = 'results/bbs-pred-results.R')
