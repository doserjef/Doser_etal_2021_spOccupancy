# spIntPGOcc-sim.R: this script simulates data for use in a spatially-explicit
#                   integrated occupancy model across 40,000 locations. We fit
#                   the model using the spIntPGOcc function in spOccupancy. The 
#                   first data source has low detection but is widespread, 
#                   second data source has medium detection and is somewhat 
#                   widespread, third data source is of the highest quality but 
#                   is spatially limited. 

# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(coda)
library(spOccupancy)
library(tidyverse)
library(viridis)
library(ggpubr)
library(ggthemes)
library(scales)
set.seed(10101)

# Simulate Data -----------------------------------------------------------
# Number of locations in each direction. This is the total region of interest
# where some sites may or may not have a data source. 
J.x <- 200
J.y <- 200
J.all <- J.x * J.y
# Number of data sources.
n.data <- 3
# Sites for each data source. 
J.obs <- c(25000, 15000, 5000)
# Replicates for each data source.
n.rep <- list()
n.rep[[1]] <- sample(1:2, size = J.obs[1], replace = TRUE)
n.rep[[2]] <- sample(2:3, size = J.obs[2], replace = TRUE)
n.rep[[3]] <- sample(2:4, size = J.obs[3], replace = TRUE)
# Occupancy covariates
beta <- c(0, -0.5, 1.0)
# Number of occurrence covariates
p.occ <- length(beta)
# Detection covariates
alpha <- list()
alpha[[1]] <- c(-1, 0.4)
alpha[[2]] <- c(0, -0.5)
alpha[[3]] <- c(1, 0.8)
# Number of detection covariates in each data source
p.det.long <- sapply(alpha, length)
# Total number of detection covariates
p.det <- sum(p.det.long)
# Spatial parameters
sigma.sq <- 2
phi <- 3 / .6
sp <- TRUE

# This is currently commented out as it takes a few hours to simulate the 
# data set as the spatial random effects are drawn from a full Gaussian process. 
# The simulated data set is read in below the commented code.
# Simulate occupancy data.
# dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs, 
# 		 n.rep = n.rep, beta = beta, alpha = alpha, sp = sp, cov.model = 'exponential', 
# 		 sigma.sq = sigma.sq, phi = phi)
# Save and load -----------------------------------------------------------
#save(data, file = 'data/simulation-data.rda')
load("simulation-data.rda")

# Package data for spIntPGOcc ---------------------------------------------
y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites
X.0 <- dat$X.pred
psi.0 <- dat$psi.pred
coords <- as.matrix(dat$coords.obs)
coords.0 <- as.matrix(dat$coords.pred)

# Package all data into a list
occ.covs <- X[, -1, drop = FALSE]
colnames(occ.covs) <- c('occ.cov.1', 'occ.cov.2')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2])
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2]) 
det.covs[[3]] <- list(det.cov.3.1 = X.p[[3]][, , 2])
data.list <- list(y = y, 
		  occ.covs = occ.covs,
		  det.covs = det.covs, 
		  sites = sites, 
		  coords = coords)

J <- length(dat$z.obs)

inits.list <- list(alpha = list(0, 0, 0), 
		      beta = c(0, -0.3, 0.7), 
		      phi = 3 / .5, 
		      sigma.sq = 2, 
		      w = rep(0, J), 
		      z = rep(1, J))

prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
		   alpha.normal = list(mean = list(0, 0, 0), 
			               var = list(2.72, 2.72, 2.72)),
		   phi.unif = c(3/1, 3/.1), 
		   sigma.sq.ig = c(2, 2)) 
# Tuning
tuning.list <- list(phi = 0.1) 

# Number of batches
n.batch <- 2000
# Batch length
batch.length <- 25

# Fit the model
out <- spIntPGOcc(occ.formula = ~ occ.cov.1 + occ.cov.2, 
		  det.formula = list(f.1 = ~ det.cov.1.1, 
				     f.2 = ~ det.cov.2.1, 
				     f.3 = ~ det.cov.3.1), 
		  data = data.list,  
		  inits = inits.list, 
		  n.batch = n.batch, 
		  batch.length = batch.length, 
		  accept.rate = 0.43, 
		  priors = prior.list, 
		  cov.model = "exponential", 
		  tuning = tuning.list, 
		  n.omp.threads = 10, 
		  verbose = TRUE, 
		  NNGP = TRUE, 
		  n.neighbors = 5, 
		  search.type = 'cb', 
		  n.report = 120, 
		  n.burn = 20000, 
		  n.thin = 30) 
# Run time
out$run.time # Approximately 78 min
save(out, file = 'results/sim-samples.rda')

# Summary stats
summary(out)

# Predict at non-sampled locations
out.pred <- predict(out, X.0, coords.0, n.omp.threads = 10)

# Plot of fitted/predicted vs. true occurrence values ---------------------
psi.quants <- apply(out$psi.samples, 2, quantile, c(0.025, 0.5, 0.975))
psi.pred.quants <- apply(out.pred$psi.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
psi.true <- dat$psi.obs
psi.0.true <- dat$psi.pred
coords.all <- rbind(coords, coords.0)
colnames(coords.all) <- c("Easting", "Northing")
plot.df <- data.frame(coords.all, 
		      psi.true = c(psi.true, psi.0.true), 
		      psi.fit = c(psi.quants[2, ], psi.pred.quants[2, ]))
psi.quants <- apply(out$psi.samples, 2, quantile, c(0.025, 0.5, 0.975))
psi.sd <- apply(out$psi.samples, 2, sd)
psi.pred.quants <- apply(out.pred$psi.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
psi.pred.sd <- apply(out.pred$psi.0.samples, 2, sd)
psi.true <- dat$psi.obs
psi.0.true <- dat$psi.pred
w.quants <- apply(out$w.samples, 2, quantile, c(0.025, 0.5, 0.975))
w.pred.quants <- apply(out.pred$w.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
w.true <- dat$w.obs
w.0.true <- dat$w.pred
w.sd <- apply(out$w.samples, 2, sd)
w.pred.sd <- apply(out.pred$w.0.samples, 2, sd)
coords.all <- rbind(coords, coords.0)
colnames(coords.all) <- c("Easting", "Northing")
plot.df <- data.frame(coords.all, 
		      psi.true = c(psi.true, psi.0.true), 
		      psi.fit = c(psi.quants[2, ], psi.pred.quants[2, ]), 
		      psi.sd = c(psi.sd, psi.pred.sd), 
		      w.true = c(w.true, w.0.true), 
		      w.fit = c(w.quants[2, ], w.pred.quants[2, ]), 
		      w.sd = c(w.sd, w.pred.sd))
true.psi.plot <- ggplot(data = plot.df, aes(x = Easting, y = Northing, fill = psi.true)) + 
  geom_raster() + 
  scale_fill_viridis_c(labels = label_number(accuracy = 0.01)) + 
  #scale_fill_gradientn(colors = magma(10)) + 
  theme_bw(base_size = 18) +
  guides(fill = FALSE) + 
  labs(title = "(A) True", fill = expression(psi))

fit.psi.plot <- ggplot(data = plot.df, aes(x = Easting, y = Northing, fill = psi.fit)) + 
  geom_raster() + 
  scale_fill_viridis_c(labels = label_number(accuracy = 0.01)) + 
  theme_bw(base_size = 18) + 
  labs(title = "(B) Predicted", fill = expression(psi)) 

sd.psi.plot <- ggplot(data = plot.df, aes(x = Easting, y = Northing, fill = psi.sd)) + 
  geom_raster() + 
  scale_fill_viridis_c(labels = label_number(accuracy = 0.01)) + 
  theme_bw(base_size = 18) + 
  labs(title = "(C) Standard Deviation", fill = 'SD') 

true.w.plot <- ggplot(data = plot.df, aes(x = Easting, y = Northing, fill = w.true)) + 
  geom_raster() + 
  scale_fill_gradientn(colors = magma(10), labels = label_number(accuracy = 0.01)) + 
  theme_bw(base_size = 18) +
  guides(fill = FALSE) + 
  labs(title = "(C) True", fill = "w")

fit.w.plot <- ggplot(data = plot.df, aes(x = Easting, y = Northing, fill = w.fit)) + 
  geom_raster() + 
  scale_fill_gradientn(colors = magma(10), labels = label_number(accuracy = 0.1)) + 
  theme_bw(base_size = 18) + 
  labs(title = "(D) Predicted", fill = "w") 

full.plot <- ggarrange(true.psi.plot, fit.psi.plot, true.w.plot, fit.w.plot, nrow = 2, ncol = 2, 
		       widths = c(4.5, 5.5)) 
# Save image if desired
# ggsave(filename = "sim-psi.pdf", plot = full.plot, height = 12, width = 12, units = 'in')

# Standard deviations
sd.psi.plot <- ggplot(data = plot.df, aes(x = Easting, y = Northing, fill = psi.sd)) + 
  geom_raster() + 
  scale_fill_viridis_c(labels = label_number(accuracy = 0.01)) + 
  theme_bw(base_size = 18) + 
  labs(title = "(A) Occurrence Probability", fill = 'SD') 

sd.w.plot <- ggplot(data = plot.df, aes(x = Easting, y = Northing, fill = w.sd)) + 
  geom_raster() + 
  scale_fill_gradientn(colors = magma(10), labels = label_number(accuracy = 0.1)) + 
  theme_bw(base_size = 18) + 
  labs(title = "(B) Spatial Process", fill = 'SD') 

# Not in the manuscript
full.sd.plot <- ggarrange(sd.psi.plot, sd.w.plot, nrow = 1, ncol = 2) 
ggsave(filename = "sim-psi-sd.pdf", plot = full.sd.plot, height = 6, width = 12, units = 'in')

# Plot of the Data Point Locations ----------------------------------------
# Not in the manuscript
plot.df$ds.1 <- 0
plot.df$ds.2 <- 0
plot.df$ds.3 <- 0
plot.df[sites[[1]], 'ds.1'] <- 1
plot.df[sites[[2]], 'ds.2'] <- 1
plot.df[sites[[3]], 'ds.3'] <- 1
plot.df$ds.locs <- 0
plot.df$ds.locs <- ifelse(plot.df$ds.1 == 1 & plot.df$ds.2 == 0 & plot.df$ds.3 == 0, 
			  1, plot.df$ds.locs)
plot.df$ds.locs <- ifelse(plot.df$ds.1 == 0 & plot.df$ds.2 == 1 & plot.df$ds.3 == 0, 
			  1, plot.df$ds.locs)
plot.df$ds.locs <- ifelse(plot.df$ds.1 == 0 & plot.df$ds.2 == 0 & plot.df$ds.3 == 1, 
			  1, plot.df$ds.locs)
plot.df$ds.locs <- ifelse(plot.df$ds.1 == 1 & plot.df$ds.2 == 1 & plot.df$ds.3 == 0, 
			  2, plot.df$ds.locs)
plot.df$ds.locs <- ifelse(plot.df$ds.1 == 1 & plot.df$ds.2 == 0 & plot.df$ds.3 == 1, 
			  2, plot.df$ds.locs)
plot.df$ds.locs <- ifelse(plot.df$ds.1 == 0 & plot.df$ds.2 == 1 & plot.df$ds.3 == 1, 
			  2, plot.df$ds.locs)
plot.df$ds.locs <- ifelse(plot.df$ds.1 == 1 & plot.df$ds.2 == 1 & plot.df$ds.3 == 1, 
			  3, plot.df$ds.locs)
sampled.plot.df <- plot.df %>%
  filter(ds.locs != 0)
sampled.plot.df$ds.locs <- factor(sampled.plot.df$ds.locs, levels = c(1, 2, 3), 
			          order = TRUE)

points.plot <- ggplot(data = sampled.plot.df, aes(x = Easting, y = Northing, fill = ds.locs)) + 
  geom_raster() + 
  scale_fill_colorblind() + 
  theme_bw(base_size = 18) + 
  labs(fill = 'Number of \nData Sources')
ggsave(filename = "sim-points-location.pdf", plot = points.plot, height = 6, width = 8, units = 'in')
