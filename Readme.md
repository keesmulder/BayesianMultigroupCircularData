# Bayesian Multigroup Circular Data
## Kees Mulder, 2014

This project contains all code for the paper "Extending Bayesian analysis of circular data to multiple groups", Mulder, K.T., & Klugkist, I. (2014). It may be used to replicate all tables and data contained within.  


R Package dependencies:

- Rcpp
- BH
- xtable
- ggplot2
- abind

Additionally, the C++ library 'boost' must be installed on the system, as well as Rtools.

## Usage

The simulation is performed from Simulation/runSimulationStudy.R, which calls all more elementary functions. For application purposes, the MCMC samplers are found in DataAnalysis/(Gibbs, MH, Rejection), and can be used by sourcing DW.R, VMMH.R and FM.R, respectively. 

# Files

## Data

Generate data that will be analyzed. Datasets will be saved in the folder Data/Datasets/.  To save space, actual data is not added to this archive. 

### `generateData.R` 

Generates a number of datasets with given properties and places them in the Datasets folder. 

### `rvmc.cpp`

An Rcpp function used to generate a single sample of von Mises distributed data quickly. 

## DataAnalysis

### `DW.R, VMMH.R, FM.R`

Three different MCMC methods are provided in the following folders: Gibbs, MH, and Rejection. Each has an R function to set up the sampler, and an Rcpp method for the core method. 

### `describeCirc.R`

Contains basic functions to analyze circular data. For more information, see annotation within this file. 

### `VenterMode.cpp`

Rcpp algorithm to calculate the highest density interval (`hmodeci()`) and the mode, which is the mean of  the (`hmode()`). This uses some bandwith, here called `cip`. This is according to Venter (1967).


## Simulation

Contains all the code necessary to analyze the generated data repeatedly with some MCMC sampler. 

### `SimulateVM.R`

Simulates a single cell.

### `SimulationStudyVM.R`

Runs SimulateVM.R for each cell.

### `runSimulationStudy.R`

Runs Simulation study for each of the different methods, and for either 1 or three groups. 

### `vonMises.R`

Contains some methods that describe the von Mises distribution. In particular, an approximation of kappa given in Fisher (1995) is used. 

## SimulationAnalysis

### `AnalysisHelperFunctions.R` 

Contains functions that are used in generating tables.

### `GenerateTables.R` 
Generates the tables as in the paper. 

### `ExampleRunsDWMHFM.R, IllustrateRandMu.R` 
Generate the figures as in the paper.

## Spread
Contains ways the research was spread, presented: articles and slides, as well as figures and tables. 

## Supplement
Contains files that explain the inner workings of the samplers, and some derivations regarding them. 

 
