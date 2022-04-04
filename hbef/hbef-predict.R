# hbef-predict.R: this script predicts occurrence for all twelve 
#                 foliage-gleaning birds across the Hubbard Brook Experimental
#                 Forest. The latent occurrence is then summarized to predict
#                 species richness across the forest. 
rm(list = ls())
library(spOccupancy)
library(coda)

# Read in non-spatial model 
load("results/hbef-msPGOcc-3-2500-2022-03-21.R")

elev <- (hbefElev$val - mean(hbef2015$occ.covs[, 1])) / sd(hbef2015$occ.covs[, 1])
X.0 <- cbind(1, elev, elev^2)
out.pred <- predict(out, X.0)

#save(out.pred, 'results/hbef-pred-msPGOcc.R')
rich.samples <- apply(out.pred$z.0.samples, c(1, 3), sum)
rich.mean <- apply(rich.samples, 2, mean)
rich.sd <- apply(rich.samples, 2, sd)

save(rich.mean, rich.sd, file = 'results/hbef-msPGOcc-pred-richness.R')
