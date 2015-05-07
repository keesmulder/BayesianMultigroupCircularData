# ----------------------------------------------------------
#   DW_chosen_k_sampler.R
# MCMC-Sampling algorithm for circular data, following Damien & Walker (1999).
# This is a specified version for testing which k is retained.
#
# Kees Tim Mulder
# Last updated: November 2014
#
# This work was supported by a Vidi grant awarded to I. Klugkist from the
# Dutch Organization for Scientific research (NWO 452-12-010).
# ----------------------------------------------------------


require(Rcpp)
require(BH)
Rcpp::sourceCpp('Simulation/DW_chosen_k_sampler.cpp')
Rcpp::sourceCpp('Data/rvmc.cpp')

# FUNCTION DEFINITIONS ------------------------------------------------------
exp_cum          <- function(x, rate = 1)     1 - exp(- rate * x)
exp_inv_cum      <- function(x, rate = 1)     -(log(1-x)/rate)
getLambda        <- function (k) factorial(k)^-2 * 0.5^(2*k)


# SAMPLER ----------------------------------------------------------------
# Burn also gets multiplied by lag.
DW <- function(th, R_0=rep(0, length(th)), mu_0=rep(0, length(th)),
                         c=rep(0, length(th)), Q=10000, burn=0, lag = 1,
                         mu_start = 0, kp_start = 2, start_w = 4,
                         Z = 85) {
  if(!is.list(th)) stop("This function requires the data to be entered as a list of groups of angles.")

  # print(th)

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
  lambda    <- getLambda(1:Z)

  # Run the Gibbs sampler.
  out <- DWC(mu_start, start_w, kp_start,
                            lambda, mu_n, R, R_t, m,
                            Z, Qb+1, lag)

  chosenk <- out[, 1]
  kappa <- out[-(1:burn), 2]
  mu    <- as.matrix(out[-(1:burn), 3:ncol(out)])%% (2*pi)


  # Return the results.
  list(mu = mu,
       kappa = kappa,
       spec = list(mu_n = mu_n,
                   R_t = R_t,
                   chosenk = chosenk)
  )
}








