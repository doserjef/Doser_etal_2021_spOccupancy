# hbef-msPGOcc.R: this script summarizes results from the Hubbard
#                 Brook Experimental Forest case study. The script produces
#                 tables and figures included in the manuscript.
# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(spOccupancy)
library(coda)
library(sf)
library(stars)
library(tidyverse)
library(viridis)
library(ggpubr)
library(RColorBrewer)

# Subroutines -------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}


# Read in results files ---------------------------------------------------
# NOTE: results files are too large to place on GitHub, so please contact
#       the first author (Jeff Doser doserjef@msu.edu) if you want them.
#       Also note they don't take that long to run so you can just run them
#       using the scripts in the hbef/ directory and produce results scripts
#       on your own.
# Load nonspatial MSOMs
load("results/hbef-msPGOcc-1-2500-2022-03-21.R")
out.1 <- out
load("results/hbef-msPGOcc-2-2500-2022-03-21.R")
out.2 <- out
load("results/hbef-msPGOcc-3-2500-2022-03-21.R")
out.3 <- out
# Load spatial MSOMs
load("results/hbef-spMsPGOcc-1-2022-03-21.R")
out.sp.1 <- out
load("results/hbef-spMsPGOcc-2-2022-03-21.R")
out.sp.2 <- out
load("results/hbef-spMsPGOcc-3-2022-03-21.R")
out.sp.3 <- out
# Load nonspatial MSOM intercept only 
# Object loaded is out.int
load("results/hbef-msPGOcc-int-no-cross-2022-03-21.R")
# Load spatial MSOM intercept only 
# Ojbect loaded is out.sp.int
load("results/hbef-spMsPGOcc-int-no-cross-2022-03-21.R")

# Convergence Assessment --------------------------------------------------
# Non-spatial model -------------------
# Community level
gelman.diag(mcmc.list(out.1$beta.comm.samples, out.2$beta.comm.samples, 
		      out.3$beta.comm.samples))
gelman.diag(mcmc.list(out.1$alpha.comm.samples, out.2$alpha.comm.samples, 
		      out.3$alpha.comm.samples))
gelman.diag(mcmc.list(out.1$tau.sq.beta.samples, out.2$tau.sq.beta.samples, 
		      out.3$tau.sq.beta.samples))
gelman.diag(mcmc.list(out.1$tau.sq.alpha.samples, out.2$tau.sq.alpha.samples, 
		      out.3$tau.sq.alpha.samples))
# Occurrence coefficients
gelman.diag(mcmc.list(out.1$beta.samples, out.2$beta.samples, 
		      out.3$beta.samples))
# Detection coefficients
gelman.diag(mcmc.list(out.1$alpha.samples, out.2$alpha.samples, 
		      out.3$alpha.samples))
# Spatial model -----------------------
gelman.diag(mcmc.list(out.sp.1$beta.comm.samples, out.sp.2$beta.comm.samples, 
		      out.sp.3$beta.comm.samples))
gelman.diag(mcmc.list(out.sp.1$alpha.comm.samples, out.sp.2$alpha.comm.samples, 
		      out.sp.3$alpha.comm.samples))
gelman.diag(mcmc.list(out.sp.1$tau.sq.beta.samples, out.sp.2$tau.sq.beta.samples, 
		      out.sp.3$tau.sq.beta.samples))
gelman.diag(mcmc.list(out.sp.1$tau.sq.alpha.samples, out.sp.2$tau.sq.alpha.samples, 
		      out.sp.3$tau.sq.alpha.samples))
# Occurrence coefficients
gelman.diag(mcmc.list(out.sp.1$beta.samples, out.sp.2$beta.samples, 
		      out.sp.3$beta.samples))
# Detection coefficients
gelman.diag(mcmc.list(out.sp.1$alpha.samples, out.sp.2$alpha.samples, 
		      out.sp.3$alpha.samples))
# Spatial parameters
gelman.diag(mcmc.list(out.sp.1$theta.samples, out.sp.2$theta.samples, 
		      out.sp.3$theta.samples))

# Posterior Summaries -----------------------------------------------------
out <- out.3
out.sp <- out.sp.3

# Community level
summary(out, level = 'community')
summary(out.sp, level = 'community')
summary(out.int, level = 'community')
summary(out.sp.int, level = 'community')

# Species level
summary(out, level = 'species')
summary(out.sp, level = 'species')
summary(out.int, level = 'species')
summary(out.sp.int, level = 'species')

# Posterior predictive checks ---------------------------------------------
ppc.out.sites <- ppcOcc(out, 'freeman-tukey', group = 1)
ppc.out.rep <- ppcOcc(out, 'freeman-tukey', group = 2)
ppc.out.sp.sites <- ppcOcc(out.sp, 'freeman-tukey', group = 1)
ppc.out.sp.rep <- ppcOcc(out.sp, 'freeman-tukey', group = 2)
# BPVs across replicates --------------
summary(ppc.out.rep, level = 'both')
summary(ppc.out.sp.rep, level = 'both')
# BPVs across sites -------------------
summary(ppc.out.sites, level = 'both')
summary(ppc.out.sp.sites, level = 'both')

# Model comparison --------------------------------------------------------
waicOcc(out)
waicOcc(out.sp)
waicOcc(out.int)
waicOcc(out.sp.int)

# Compare cross-validation results
load("results/hbef-msPGOcc-cross-val-3-1000-2022-03-25.R")
out.ms.cross <- out
load("results/hbef-spMsPGOcc-cross-val-3-2022-03-25.R")
out.sp.ms.cross <- out
# Intercept only model
load("results/hbef-msPGOcc-int-2022-03-21.R")
# Spatial intercept only model
load("results/hbef-spMsPGOcc-int-2022-03-25.R")
sum(out.ms.cross$k.fold.deviance)
sum(out.sp.ms.cross$k.fold.deviance)
sum(out.int$k.fold.deviance)
sum(out.sp.int$k.fold.deviance)

# Read in predictions -----------------------------------------------------
# Only looking at non-spatial predictions since they are better according to 
# k-fold cross-validation. 
load("results/hbef-msPGOcc-pred-richness.R")

rich.df <- data.frame(Easting = hbefElev$Easting, 
		      Northing = hbefElev$Northing, 
		      elev = hbefElev$val,
		      rich.mean = rich.mean, 
		      rich.sd = rich.sd)

coords.sf <- st_as_sf(x = rich.df, 
		      coords = c('Easting', 'Northing'), 
		      crs = "+proj=utm +zone=19 +units=m +datum=NAD83")

hbef <- st_read('hbef/hbef-spatial/', 'hbef_wsheds')
# Join all watersheds together
hbef <- st_union(hbef)
coords.stars <- st_as_stars(rich.df)

# Mean richness
rich.mean.plot <- ggplot(data = hbef) + 
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = rich.mean)) + 
  geom_sf(alpha = 0, col = 'gray58') + 
  scale_fill_viridis_c(na.value = 'transparent') +
  theme_bw(base_size = 14) + 
  labs(fill = 'Mean Species\nRichness', x = 'Longitude', y = 'Latitude') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Elevation
elev.plot <- ggplot(data = hbef) + 
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = elev)) + 
  geom_sf(alpha = 0, col = 'gray58') + 
  scale_fill_viridis_c(na.value = 'transparent') +
  theme_bw(base_size = 14) + 
  labs(fill = 'Elevation', x = 'Longitude', y = 'Latitude') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../figures/FigS1.pdf", height = 3.5, width = 7, units = 'in')

# Standard deviation richness
rich.sd.plot <- ggplot(data = hbef) + 
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = rich.sd)) + 
  geom_sf(alpha = 0, col = 'gray58') + 
  scale_fill_viridis_c(na.value = 'transparent') +
  theme_bw(base_size = 14) + 
  labs(fill = 'SD Species\nRichness', x = 'Longitude', y = 'Latitude') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Figure 2
ggarrange(rich.mean.plot, rich.sd.plot, ncol = 1, nrow = 2, labels = c('(A)', '(B)'))
ggsave("../figures/Fig2.pdf", height = 7, width = 7, units = 'in')


# Comparison across different numbers of neighbors ------------------------
# Read in all spatial MSOMs fit with different number of neighbors, which is
# denoted by the number in the saved object. 
load("results/hbef-spMsPGOcc-nn-1-1-2022-03-25.R")
out.1 <- out
load("results/hbef-spMsPGOcc-nn-3-3-2022-03-25.R")
out.3 <- out
load("results/hbef-spMsPGOcc-nn-5-3-2022-03-25.R")
out.5 <- out
load("results/hbef-spMsPGOcc-nn-7-3-2022-03-25.R")
out.7 <- out
load("results/hbef-spMsPGOcc-nn-9-3-2022-03-25.R")
out.9 <- out
load("results/hbef-spMsPGOcc-nn-11-3-2022-03-25.R")
out.11 <- out
load("results/hbef-spMsPGOcc-nn-13-3-2022-03-25.R")
out.13 <- out
load("results/hbef-spMsPGOcc-nn-15-3-2022-03-25.R")
out.15 <- out
load("results/hbef-spMsPGOcc-nn-17-3-2022-03-25.R")
out.17 <- out
load("results/hbef-spMsPGOcc-nn-19-3-2022-03-25.R")
out.19 <- out
load("results/hbef-spMsPGOcc-nn-21-3-2022-03-25.R")
out.21 <- out

# Put them all in a list
out.list <- list(out.1, out.3, out.5, out.7, out.9, out.11, out.13, out.15, out.17, 
		 out.19, out.21)
run.time.all <- lapply(out.list, function(a) a$run.time)
waic.out.all <- lapply(out.list, waicOcc)
beta.comm.mean.all <- lapply(out.list, function(a) apply(a$beta.comm.samples, 2, mean))
beta.comm.low.all <- lapply(out.list, function(a) apply(a$beta.comm.samples, 2, quantile, 0.025))
beta.comm.high.all <- lapply(out.list, function(a) apply(a$beta.comm.samples, 2, quantile, 0.975))
# Plot all Results --------------------------------------------------------
plot.dat <- data.frame(run.time = sapply(run.time.all, function(a) a[3]) / 60, 
		       n.neighbors = sapply(out.list, function(a) a$n.neighbors), 
		       waic = sapply(waic.out.all, function(a) a[3]), 
		       beta.0.mean = sapply(beta.comm.mean.all, function(a) a[1]),
		       beta.1.mean = sapply(beta.comm.mean.all, function(a) a[2]),
		       beta.2.mean = sapply(beta.comm.mean.all, function(a) a[3]))

# Plot of run time (Figure 3a)
run.time.plot <- ggplot(data = plot.dat, aes(x = n.neighbors)) + 
  geom_line(aes(y = run.time), lty = 3) + 
  geom_point(aes(y = run.time), fill = 'lightskyblue1', pch = 21, size = 4) + 
  theme_bw(base_size = 18) + 
  labs(x = 'Number of Neighbors', y = 'Run Time (min)')

new.plot.dat <- data.frame(n.neighbors = rep(plot.dat$n.neighbors, 3), 
			   beta.mean = c(plot.dat$beta.0.mean, 
					 plot.dat$beta.1.mean, 
					 plot.dat$beta.2.mean), 
			   param = rep(c('Intercept', 'Elevation', 'Elevation.2'), 
				       each = nrow(plot.dat)))
# Plot of community-level estimates (Figure 3b)
coefs.plot <- ggplot(data = new.plot.dat, aes(x = n.neighbors, y = beta.mean, fill = param)) + 
  geom_point(pch = 21, size = 4) + 
  scale_fill_viridis_d(labels = c('Intercept', 'Elevation', expression(Elevation^2))) + 
  theme_bw(base_size = 18) + 
  theme(legend.text.align = 0, 
	legend.position = c(0.79, 0.3), 
	legend.background = element_rect(fill = 'white', color = 'gray')) + 
  labs(x = 'Number of Neighbors', y = 'Estimate', fill = 'Parameter')

# Figure S2
full.plot <- ggarrange(run.time.plot, coefs.plot, nrow = 1, ncol = 2, labels = c('(A)', '(B)'))
ggsave('../figures/hbef-nn-comparison.pdf', height = 5, width = 10, units = 'in')
