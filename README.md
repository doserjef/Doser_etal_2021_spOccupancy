# [spOccupancy: An R package for single species, multispecies, and integrated spatial occupancy models](https://arxiv.org/pdf/2111.12163.pdf)

### In Review

### Jeffrey W. Doser, Andrew O. Finley, Marc K&eacute;ry, Elise F. Zipkin 

### Package [Website](https://www.jeffdoser.com/files/spoccupancy-web/) and [Repository](https://github.com/doserjef/spOccupancy/)

### Please contact the first author for questions about the code or data used in the empirical case studies: Jeffrey W. Doser (doserjef@msu.edu)

---------------------------------

## Abstract

1. Occupancy modeling is a common approach to estimate spatial and temporal species distribution patterns, while explicitly accounting for measurement errors common in detection-nondetection data. Numerous extensions of the basic single species occupancy model exist to address dynamics, multiple species or states, interactions, false positive errors, autocorrelation, and to integrate multiple data sources, yet specialized and computationally efficient software to fit spatial models to large data sets is scarce or absent. 
2. We introduce the `spOccupancy` `R` package designed to fit single species, multispecies, and integrated spatially-explicit occupancy models. Using a Bayesian framework, we leverage P&oacute;ly Gamma data augmentation and Nearest Neighbor Gaussian Processes to ensure models are computationally efficient for potentially massive data sets.
3. `spOccupancy` provides computationally efficient and user-friendly functions for data simulation, model fitting, model validation (by posterior predictive checks), model comparison (using information criteria and k-fold cross-validation), and out-of-sample prediction. We illustrate the package's functionality via a vignette, simulated data analysis, and two bird case studies, in which we estimate occurrence of the Black-throated Green Warbler (*Setophaga virens*) across the eastern USA and species richness of a foliage-gleaning bird community in the Hubbard Brook Experimental Forest in New Hampshire, USA. 
4. The `spOccupancy` package provides a user-friendly approach to fit a variety of single and multispecies occupancy models, making it straightforward to address detection biases and spatial autocorrelation in these specialized species distribution models even for large data sets.  

## Repository Directory

All code and resulting model objects were created and saved using spOccupancy v0.3.0.

### [BBS](./bbs)

Contains all code and data for case study of the Black-throated Green Warbler distribution across the eastern USA. 

+ `data`: directory containing the raw BBS data used in the analysis.
+ `bbs-PGOcc-cross-val.R`: script to run nonspatial single species occupancy model with cross-validation.
+ `bbs-PGOcc.R`: script to run nonspatial single species occupancy model.
+ `bbs-data-prep.R`: script to prepare the raw data in the `data` subdirectory for analysis in `spOccupancy`. 
+ `bbs-pred-data-prep.R`: script to prepare the raw data in the `data` subdirectory for prediction across the eastern US in `spOccupancy`.
+ `bbs-predict.R`: code to predict occurrence across the eastern US. 
+ `bbs-spPGOcc-GP-cross-val.R`: script to run spatial single species occupancy model using a full Gaussian process with cross-validation.
+ `bbs-spPGOcc-cross-val.R`: script to run spatial single species occupancy model using an NNGP with cross-validation.
+ `bbs-spPGOcc.R`: script to run spatial single species occupancy model using an NNGP.
+ `bbs-spPGOccGP.R`: script to run spatial single species occupancy model using a full Gaussian process.
+ `pfile-1`, `pfile-2`, `pfile-3`, `pfile-sp-1`, `pfile-sp-2`, `pfile-sp-3`: files used to specify initial values when running files across multiple cores from the command line.
+ `summary-bbs.R`: script to perform summary analyses of all model results. Code to produce Figure 1 and Table 2 is in this script.

### [HBEF](./hbef)

Contains all code and data for case study of the foliage-gleaning bird community in the Hubbard Brook Experimental Forest.

+ `hbef-spatial`: directory containing shapefiles for creation of Figure 2 in the manuscript.
+ `hbef-msPGOcc-cross-val.R`: script to run nonspatial multispecies occupancy model with cross-validation.
+ `hbef-msPGOcc-int.R`: script to run nonspatial multispecies occupancy model with an intercept only model for occurrence.
+ `hbef-msPGOcc.R`: script to run nonspatial multispecies occupancy model.
+ `hbef-predict.R`: code to predict species-specific occurrence and species richness from a nonspatial multispecies occupancy model.
+ `hbef-spMsPGOcc-cross-val.R`: script to run spatial multispecies occupancy model with cross-validation.
+ `hbef-spMsPGOcc-int.R`: script to run spatial multispecies occupancy model with an intercept only model for occurrence.
+ `hbef-spMsPGOcc-nn.R`: script to compare spatial multispecies occupancy models fit with different numbers of nearest neighbors.
+ `hbef-spMsPGOcc.R`: script to run spatial multispecies occupancy model.
+ `pfile-1`, `pfile-2`, `pfile-3`, `pfile-sp-1`, `pfile-sp-2`, `pfile-sp-3`: files used to specify initial values when running files across multiple cores from the command line.
+ `summary-hbef.R`: script to perform summary analyses of all model results. Code to produce Figure 2, Table 3, and Figure 3 is in this script.

### [Simulations](./simulations)

Contains code and data for analysis of a simulated data set. 

+ `simulation-data.rda`: simulated data set obtained using `simIntOcc`. 
+ `spIntPGOcc-sim.R`: script to fit spatial integrated occupancy model for the simulated data set and analyze the results. Code to produce Figure 4 is in this script.




