(1) Contains the main article. 

(2) Contains code to generate the datasets as used in the paper. They will be saved in the folder.  To save space, actual data is not added to this archive. 

(3) Three different MCMC methods are provided. They require Rcpp.

(4) Code for simulation is provided. 
-- SimulateVM.R simulates a single cell.
-- SimulationStudyVM.R runs SimulateVM.R for each cell.
-- runSimulationStudy runs Simulation study for each of the different methods, and for either 1 or three groups. 
-- AnalysisHelperFunctions.R provides methods to extract information from the simulation results.
-- VenterMode.cpp is a short RCPP programm to compute the mode according to Venter (1967)

(5) Code for analysis
-- AnalysisHelperFunctions.R contains functions that are used in generating tables.
-- GenerateTables generates the tables as in the paper. 
-- ExampleRunsDWMHFM and IllustrateRandMu generate the figures as in the paper.

(6) Basic Codes
-- DescribeCirc.R contains a collection of functions that describe or plot circular datasets. For more information, see annotations within that file. 
-- rvmc.cpp is a generally usable code that will generate von Mises random variates.


(7) Supplemental Files
-- Contains files that explain the inner workings of the samplers, and some derivations regarding them. 

 