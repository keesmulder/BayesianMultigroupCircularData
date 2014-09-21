source("(4) Code for simulation/SimulateVM.R")

# Perform a simulation study on a sampler with given properties.
simulationStudyVM <- function (samplername, nsim, ns=c(5, 30, 100), kappas=c(0.1, 1, 4, 16, 32), 
                               meandifs=20*(pi/180), J=1, 
                               Q=10000, burn=500, lag=1, 
                               printcell=TRUE, printsim=FALSE, ...) {
  
  
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
  # for n, we want to know if it breaks down.
  # for kappa, we want to see what the limits are
  FSR <- array(NA, dim = fd, 
               dimnames = list(simnr = 1:nsim, 
                               outcomes=1:noutcome,
                               n = ns, 
                               kappa=kappas, 
                               meandif=meandifs))
  
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
  
  ssres <- list(Design=list(),
                Results=FSR)
  
  timeFinished <- Sys.time()
  timeDif <- timeFinished - timeStarted
  
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
  
  filename <- paste("[SimResultVM_v6", paste0("nsim", nsim), samplername, 
                    paste0("n", paste(ns, collapse=",")), 
                    paste0("k", paste(kappas, collapse=",")), 
                    paste0("mudif", paste(round(meandifs,2), collapse=",")), 
                    paste0("J", J), paste0("Q", Q), paste0(textTime,"]"), sep = "]__[")
  
  save(ssres, file=paste0("(4) Code for simulation/Results", filename, ".rda"))
  
  
  return(ssres)
}


