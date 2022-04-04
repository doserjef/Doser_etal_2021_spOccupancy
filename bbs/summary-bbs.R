# summary-bbs.R: this script loads in result files from spOccupancy analyses
#                and compares results across models. Summary statistics and 
#                figures presented in the manuscript for the BBS case study
#                are created in this script. 
# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(spOccupancy)
library(coda)
library(tidyverse)
library(sf)
library(ggpubr)
library(stars)
library(viridis)

# Read in results files ---------------------------------------------------
# NOTE: these files are too large to place on GitHub, so please contact
#       the first author (Jeff Doser doserjef@msu.edu) if you want them. 
#       Also note they don't take that long to run so you can just run them 
#       using the scripts in the bbs/ directory and produce results scripts
#       on your own.
# Nonspatial SSOMs
load("results/bbs-PGOcc-1-2022-03-21.R")
out.1 <- out
load("results/bbs-PGOcc-2-2022-03-21.R")
out.2 <- out
load("results/bbs-PGOcc-3-2022-03-21.R")
out.3 <- out
# Spatial SSOMs with NNGP
load("results/bbs-spPGOcc-1-2022-03-22.R")
out.sp.1 <- out
load("results/bbs-spPGOcc-2-2022-03-22.R")
out.sp.2 <- out
load("results/bbs-spPGOcc-3-2022-03-22.R")
out.sp.3 <- out
# Spatial SSOMs with full Gaussian process (GP)
load("results/bbs-spPGOcc-GP-1-2022-03-22.R")
out.sp.gp.1 <- out
load("results/bbs-spPGOcc-GP-2-2022-03-22.R")
out.sp.gp.2 <- out
load("results/bbs-spPGOcc-GP-3-2022-03-22.R")
out.sp.gp.3 <- out

# Convergence Assessment --------------------------------------------------
# Nonspatial model --------------------
# Occurrence coefficients
gelman.diag(mcmc.list(out.1$beta.samples, out.2$beta.samples,
		      out.3$beta.samples))
# Detection coefficients
gelman.diag(mcmc.list(out.1$alpha.samples, out.2$alpha.samples,
		      out.3$alpha.samples))
gelman.diag(mcmc.list(out.1$sigma.sq.p.samples, out.2$sigma.sq.p.samples, 
		      out.3$sigma.sq.p.samples))
gelman.diag(mcmc.list(out.1$alpha.star.samples, out.2$alpha.star.samples, 
		      out.3$alpha.star.samples))
# NNGP Spatial model ------------------
# Occurrence coefficients
gelman.diag(mcmc.list(out.sp.1$beta.samples, out.sp.2$beta.samples,
		      out.sp.3$beta.samples))
plot(mcmc.list(out.sp.1$beta.samples, out.sp.2$beta.samples,
		      out.sp.3$beta.samples))
# Detection coefficients
gelman.diag(mcmc.list(out.sp.1$alpha.samples, out.sp.2$alpha.samples,
		      out.sp.3$alpha.samples))
gelman.diag(mcmc.list(out.sp.1$sigma.sq.p.samples, out.sp.2$sigma.sq.p.samples, 
		      out.sp.3$sigma.sq.p.samples))
gelman.diag(mcmc.list(out.sp.1$alpha.star.samples, out.sp.2$alpha.star.samples, 
		      out.sp.3$alpha.star.samples))
# Spatial parameters
gelman.diag(mcmc.list(out.sp.1$theta.samples, out.sp.2$theta.samples,
		      out.sp.3$theta.samples))
plot(mcmc.list(out.sp.1$theta.samples, out.sp.2$theta.samples,
		      out.sp.3$theta.samples))
# GP Spatial model --------------------
# Occurrence coefficients
gelman.diag(mcmc.list(out.sp.gp.1$beta.samples, out.sp.gp.2$beta.samples,
		      out.sp.gp.3$beta.samples))
plot(mcmc.list(out.sp.gp.1$beta.samples, out.sp.gp.2$beta.samples,
		      out.sp.gp.3$beta.samples))
# Detection coefficients
gelman.diag(mcmc.list(out.sp.gp.1$alpha.samples, out.sp.gp.2$alpha.samples,
		      out.sp.gp.3$alpha.samples))
gelman.diag(mcmc.list(out.sp.gp.1$sigma.sq.p.samples, out.sp.gp.2$sigma.sq.p.samples, 
		      out.sp.gp.3$sigma.sq.p.samples))
gelman.diag(mcmc.list(out.sp.gp.1$alpha.star.samples, out.sp.gp.2$alpha.star.samples, 
		      out.sp.gp.3$alpha.star.samples))
# Spatial parameters
gelman.diag(mcmc.list(out.sp.gp.1$theta.samples, out.sp.gp.2$theta.samples,
		      out.sp.gp.3$theta.samples))

# Posterior summaries -----------------------------------------------------
out <- out.3
out.sp <- out.sp.3
out.gp <- out.sp.gp.3

summary(out)
summary(out.sp)
summary(out.gp)

# Posterior predictive checks ---------------------------------------------
ppc.out.sites <- ppcOcc(out, 'freeman-tukey', group = 1)
ppc.out.sp.sites <- ppcOcc(out.sp, 'freeman-tukey', group = 1)
ppc.out.gp.sites <- ppcOcc(out.gp, 'freeman-tukey', group = 1)

# BPVs across sites --------------
summary(ppc.out.sites)
summary(ppc.out.sp.sites)
summary(ppc.out.gp.sites)

# Non-spatial
ppc.df <- data.frame(fit = ppc.out.sites$fit.y, 
		     fit.rep = ppc.out.sites$fit.y.rep, 
		     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = "Fit", xlab = "True")
lines(ppc.df$fit, ppc.df$fit, col = 'black')

# Spatial
ppc.sp.df <- data.frame(fit = ppc.out.sp.sites$fit.y, 
		     fit.rep = ppc.out.sp.sites$fit.y.rep, 
		     color = 'lightskyblue1')
ppc.sp.df$color[ppc.sp.df$fit.rep > ppc.sp.df$fit] <- 'lightsalmon'
plot(ppc.sp.df$fit, ppc.sp.df$fit.rep, bg = ppc.sp.df$color, pch = 21, 
     ylab = "Fit", xlab = "True")
lines(ppc.sp.df$fit, ppc.sp.df$fit, col = 'black')

# Model comparison --------------------------------------------------------
waic.out <- waicOcc(out)
waic.out.sp <- waicOcc(out.sp)
waic.out.gp <- waicOcc(out.gp)
waic.out
waic.out.sp
waic.out.gp

load("results/bbs-spPGOcc-cross-val-2022-03-22.R")
out.sp.cross <- out
load("results/bbs-PGOcc-cross-val-2022-03-22.R")
out.cross <- out
load("results/bbs-spPGOcc-GP-cross-val-2022-03-24.R")
out.gp.cross <- out

out.cross$k.fold.deviance
out.sp.cross$k.fold.deviance
out.gp.cross$k.fold.deviance

# Model fit plots ---------------------------------------------------------
# This plot is not in the manuscript. 
load("bbs/data/bbs-btnw-bundle.R")
# Get map of states of interest
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
btnw.states <- usa %>%
  dplyr::filter(ID %in% c('connecticut', 'delaware', 'georgia', 
		   'illinois', 'indiana', 'kentucky', 'maine', 'maryland', 
		   'massachusetts', 'michigan', 'minnesota', 'north carolina', 
		   'new hampshire', 'new jersey', 'new york', 'ohio', 
		   'pennsylvania', 'rhode island', 'south carolina', 'tennessee', 
		   'vermont', 'virginia', 'wisconsin', 'west virginia'))

coords.sf <- st_as_sf(data.frame(bbs.btnw.dat$coords, val = apply(bbs.btnw.dat$y, 1, max)),
		      coords = c("Longitude", "Latitude"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Albers equal area across contiguous US.
coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Get coordinates in Albers Equal Area
coords.albers <- st_coordinates(coords.sf.albers)
bbs.btnw.dat$coords <- coords.albers

coords.sf.albers$psi.sp.mean <- apply(out.sp$psi.samples, 2, mean)
coords.sf.albers$psi.mean <- apply(out$psi.samples, 2, mean)
coords.sf.albers$psi.gp.mean <- apply(out.gp$psi.samples, 2, mean)
coords.sf.albers$psi.sp.sd <- apply(out.sp$psi.samples, 2, sd)
coords.sf.albers$psi.sd <- apply(out$psi.samples, 2, sd)
coords.sf.albers$psi.gp.sd <- apply(out.gp$psi.samples, 2, sd)

# Map of spatial occupancy probability
psi.sp <- ggplot() + 
  geom_sf(data = btnw.states) + 
  geom_sf(data = coords.sf.albers, aes(bg = psi.sp.mean), pch = 21, size = 2.5) + 
  theme_bw(base_size = 16) + 
  scale_fill_viridis_c() + 
  labs(x = "Longitude", y = "Latitude", fill = expression(psi)) +
  theme(legend.position = c(0.13, 0.25), 
	legend.background = element_rect(fill = 'white', color = 'black'))
# Map of non-spatial occupancy probability
psi.no.sp <- ggplot() + 
  geom_sf(data = btnw.states) + 
  geom_sf(data = coords.sf.albers, aes(bg = psi.mean), pch = 21, size = 2.5) + 
  theme_bw(base_size = 16) + 
  scale_fill_viridis_c() + 
  labs(x = "Longitude", y = "Latitude", fill = expression(psi)) +
  theme(legend.position = c(0.13, 0.25), 
	legend.background = element_rect(fill = 'white', color = 'black'))
ggarrange(psi.sp, psi.no.sp, nrow = 1, ncol = 2)

# Plot of predicted values ------------------------------------------------
load("results/bbs-pred-results.R")
load("bbs/data/bbs-pred-data.rda")
psi.0.mean <- apply(psi.0.samples, 2, mean)
psi.0.sd <- apply(psi.0.samples, 2, sd)
pred.df <- st_as_sf(data.frame(coords.0, psi.0.mean, psi.0.sd),
		    coords = c('X', 'Y'),
		    crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

buffer.pred.df <- st_buffer(pred.df, dist = 0)
# Get predicted values in area of interest
pred.plot.df <- st_intersection(pred.df, st_buffer(btnw.states, 0))
coords.0.ne <- st_coordinates(pred.plot.df)
pred.plot.df <- st_drop_geometry(pred.plot.df)

plot.stars.df <- data.frame(x = coords.0.ne[, 1], y = coords.0.ne[, 2],
			    psi.0.mean = pred.plot.df$psi.0.mean,
			    psi.0.sd = pred.plot.df$psi.0.sd)
pred.stars <- st_as_stars(plot.stars.df, dims = c('x', 'y'))

# As rasters
psi.pred.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = psi.0.mean),interpolate = TRUE) +
  geom_sf(data = btnw.states, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = viridis(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme(legend.position = c(0.13, 0.25),
	legend.background = element_rect(fill = 'white', color = 'black', size = 0.75))

psi.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = psi.0.sd),interpolate = TRUE) +
  geom_sf(data = btnw.states, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = viridis(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "SD") +
  theme(legend.position = c(0.13, 0.25),
	legend.background = element_rect(fill = 'white', color = 'black', size = 0.75))

# Figure 1
ggarrange(psi.pred.plot, psi.sd.plot, nrow = 1, ncol = 2, labels = c("(A)", "(B)"),
	  font.label = list(size = 25))
# Save figure to hard drive if desired.
ggsave("../figures/Fig1.pdf", height = 7.5, width = 15, units = 'in')
