# ----------------------------------------------------------
# SimulationStudyVM.R
# Runs simulateVM.R for several different scenario's, specifically with
# different n, kappa, and meandifs.
#
# Kees Tim Mulder Last updated: November 2014
#
# This work was supported by a Vidi grant awarded to I. Klugkist from the Dutch
# Organization for Scientific research (NWO 452-12-010).
# ----------------------------------------------------------


source("Simulation/SimulateVM.R")

# Perform a simulation study on a sampler with given properties.
simulationStudyVM <- function (samplername, nsim, ns=c(5, 30, 100), kappas=c(0.1, 4, 32),
                               meandifs=20*(pi/180), J=1,
                               Q=10000, burn=500, lag=1,
                               printcell=TRUE, printsim=FALSE, returnnull=TRUE, ...) {
  # FUNCTION simulationStudyVM ----------------------------------------------
  # samplername: The MCMC-sampler to be used, passed as a string.
  # nsim: The number of datasets to read in. These must have been generated
  #       beforehand.
  # ns, kappas, meandifs: Properties of the datasets to be read in for analysis.
  #                       Respectively the sample size, concentration,
  #                       and difference between the group means.
  # J: Number of groups to be used.
  # Q: The desired number of iterations to run the chosen MCMC method.
  # burn: Number of iterations to discard as burn-in.
  # lag: Number representing a thinning factor.
  #      Only 1/lag iterations will be saved.
  # printcell: Whether the name of each scenario should be printed.
  # printsim: Whether to print the number of the dataset currently in analysis.
  # returnnull: Whether to return NULL instead of the results. If TRUE, we only
  #             have the saved results.
  # ...: Further arguments to be passed to SimulateVM().
  #
  # Returns: A list of length 2 containing the results and the design, or NULL.
  # ------------------------------------------------------------------------

  # For getting computational time.
  timeStarted <- Sys.time()

  # Get the function
  FUN <- get(samplername)

  # Do a quick testrun to get colnames, dimensions.
  testrun <- simulateVM(nsim=1, n=5, kappa=kappas[1],
                        J=J, meandif=meandifs[1],
                        Q=5, burn=1, FUN=FUN)

  # The dimensions of the summary and full results
  noutcome <- ncol(testrun)

  # Dimensions of the output.
  fd <- c(nsim, noutcome, sapply(list(ns, kappas, meandifs), length))

  # Full simulation results
  # [outcomes, nsim, n, kappa, meandif, ngroup (J)]
  FSR <- array(NA, dim = fd,
               dimnames = list(simnr = 1:nsim,
                               outcomes=1:noutcome,
                               n = ns,
                               kappa=kappas,
                               meandif=meandifs))

  # Run simulateVM for each scenario.
  for (n in ns) {
    for (kappa in kappas) {
      for (meandif in meandifs) {

        if (printcell) {
          print(paste0("n=",n, ", kappa=",kappa, ", meandif=",meandif,
                       ", ngroup=", J, ", at ", Sys.time()))
        }

        FSR[ , , as.character(n),
            as.character(kappa),
            as.character(meandif)] <-
          simulateVM(nsim=nsim, n=n, kappa=kappa, J=J, meandif=meandif,
                     Q=Q, burn=burn, lag=lag, FUN=FUN, printsim=printsim, ...)
      }
    }
  }

  dimnames(FSR)[[2]] <- colnames(testrun)

  # Add design to the results matrix.
  ssres <- list(Design=list(),
                Results=FSR)

  # Save the time simulation took.
  timeFinished <- Sys.time()
  timeDif <- timeFinished - timeStarted

  # Construct output.
  ssres$Design <- list("Sampling Method" = samplername,
                       "Number of simulations"=nsim,
                       "Sample sizes"=ns,
                       "Kappas"=kappas,
                       "Mean differences"=meandifs,
                       "Number of groups"=J,
                       "Number of Iterations"=Q,
                       "Burn-in"=burn,
                       "Lag"=lag,
                       "N_k_attempts"=ifelse(samplername=="DW", 85, NA),
                       "Time started"=timeStarted,
                       "Time completed"=timeFinished,
                       "Time Difference"=timeDif)


  textTime <- gsub(":", ".", timeFinished)

  filename <- paste("[SimResultVM_v6",
                    paste0("nsim", nsim),
                    samplername,
                    paste0("n", paste(ns, collapse=",")),
                    paste0("k", paste(kappas, collapse=",")),
                    paste0("mudif", paste(round(meandifs,2), collapse=",")),
                    paste0("J", J), paste0("Q", Q), paste0(textTime,"]"), sep = "]__[")

  # Save the results to disk as a .rda file.
  save(ssres, file=paste0("Simulation/Results/", filename, ".rda"))

  if(returnnull) {
    return(NULL)
  } else {
    return(ssres)
  }
}


