# This is an implementation of the circular data Gibbs sampler in Damien &
# Walker, 1999. A more thorough explanation of the article and this code can be
# found in the file Explanation of Damien & Walker.
require(Rcpp)

Sys.setenv("PKG_CXXFLAGS" = paste0("-I", getwd()))
sourceCpp("DataAnalysis/Gibbs/DWC.cpp")

# FUNCTION DEFINITIONS -------------------------------------------------------
# CDF and Inverse CDF (Quantile function) for the exponential distribution.
exp_cum          <- function(x, rate = 1)     1 - exp(- rate * x)
exp_inv_cum      <- function(x, rate = 1)     -(log(1-x)/rate)

# Obtains lambda_k, which is used in the Bessel function.
getLambda        <- function (k) factorial(k)^-2 * 0.5^(2*k)


DW <- function(th, R_0=rep(0, length(th)), mu_0=rep(0, length(th)),
                         c=rep(0, length(th)), Q=10000, burn=0, lag = 1,
                         mu_start = 0, kp_start = 2, start_w = 4,
                         Z = 85) {
  # SAMPLER ----------------------------------------------------------------
  # th: The circular data supplied as a list with one group per item in the list.
  # mu_0, R_0, c: Prior representing c observations in direction mu_0,
  #               with resultant length R_0.
  # Q: The desired number of iterations.
  # burn: Number of iterations to discard as burn in.
  # lag: Number representing a thinning factor.
  #      Only 1/lag iterations will be saved.
  # mu_start, kp_start, start_w: Starting values of mean, concentration, and w.
  # Z: is the chosen Z.

  if(!is.list(th)) stop(paste("This function requires the data to be entered",
                              "as a list of groups of angles.")

  # Qb is the amount of iterations plus the burn-in.
  Qb <- Q + burn

  ## Data properties --------------------
  # Number of groups
  J <- length(th)

  # Properties that are separate for each group
  n <- R <- mu_n <- numeric(J)

  for (j in 1:J) {
    jdat <- th[[j]]

    # Group sample size and sum of (co-)sines in the data
    n[j] <- length(jdat)

    # Posterior properties
    C_n     <- R_0[j] * cos(mu_0[j]) + sum(cos(jdat))
    S_n     <- R_0[j] * sin(mu_0[j]) + sum(sin(jdat))
    mu_n[j] <- atan2(S_n,  C_n) %% (2*pi)
    R[j]    <- C_n/cos(mu_n[j])
  }

  # Total sample size
  m   <- sum(n) + sum(c)

  # Sum of all resultant lengths per group
  R_t <- sum(R)

  # Sampled parameters
  mu <- matrix(NA, nrow = Q, ncol = J)
  kappa <- w <- as.numeric(rep(NA, Q))
  N_k <- numeric(Z)

  ## Initialization --------------------
  mu[1, ]   <- mu_start
  w[1]      <- start_w
  kappa[1]  <- kp_start
  lambda    <- getLambda(1:Z)

  # Run the Gibbs sampler in C++.
  out <- DWC(mu_start, start_w, kp_start,
                            lambda, mu_n, R, R_t, m,
                            Z, Qb+1, lag)

  kappa <- out[-(1:burn), 2]
  mu    <- as.matrix(out[-(1:burn), 3:ncol(out)])%% (2*pi)


  # Return the results.
  list(mu = mu,
       kappa = kappa,
       spec = list(mu_n = mu_n,
                   R_t = R_t)
  )
}
