# This is an implementation of the MH-algorithm for the von Mises distribution.
Sys.setenv("PKG_CXXFLAGS"="-Ic:/Boost/include/boost-1_55")
Rcpp::sourceCpp('(3) MCMC Samplers/MH/VMMHC.cpp')

# SAMPLER ----------------------------------------------------------------
VMMH <- function(th, R_0=rep(0, length(th)), mu_0=rep(0, length(th)), 
                       c=rep(0, length(th)), Q=10000, burn=1000, lag = 1,
                       mu_start = 0, kp_start = 2) {
  
  # Qb is the amount of iterations plus the burn-in.
  Qb <- Q + burn
  
  # Accepted counter
  accepted <- 0
  
  ## Data properties --------------------
  # Number of groups
  J <- length(th)

  # Properties that are separate for each group
  n <- R_n <- mu_n <- R_0 <- mu_0 <- numeric(J)
  
  for (j in 1:J) {
    # Group sample size and sum of (co-)sines in the data
    n[j] <-  length(th[[j]])     
    C    <- sum(cos(th[[j]]))   
    S    <- sum(sin(th[[j]]))   
    
    # Posterior properties
    C_n     <- R_0[j] * cos(mu_0[j]) + C   
    S_n     <- R_0[j] * sin(mu_0[j]) + S   
    mu_n[j] <- atan2(S_n,  C_n) %% (2*pi)  
    R_n[j]    <- C_n/cos(mu_n[j])            
  }
  
  mu <- matrix(NA, nrow = Qb, ncol = J)
  kp <- as.numeric(rep(NA, Qb))
  
  # Total sample size
  m   <- n + c
  
  # Sum of all resultant lengths per group
  R_t <- sum(R_n)
  
  # Here, we call the function that is written in C, which is in the file VMMHC.cpp.
  out      <- VMMHC(kp_start=kp_start, mu_n=mu_n, R_n=R_n, m=m, R_t=R_t, Q=Qb, lag=lag)
  
  # Rearrange output to return the required information. 
  sam      <- out$sam
  accepted <- out$acc
  
  if (burn > 0) {
    kp <- sam[-(1:burn),1]
    mu <- sam[-(1:burn),-1, drop=FALSE] %% (2*pi)  
  } else {
    kp <- sam[ , 1]
    mu <- sam[ ,-1, drop=FALSE] %% (2*pi)  
  }
  
  return(list(mu=mu,
              kappa=kp,
              spec=list(acceptance = accepted/(Qb*lag), R_t=R_t)))
}
