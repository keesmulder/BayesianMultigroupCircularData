# ----------------------------------------------------------
# SimulateVM.R
# Given some properties of data, read in all nsim datasets and analyze
# these with a given sampling function.
#
# Kees Tim Mulder
# Last updated: November 2014
#
# This work was supported by a Vidi grant awarded to I. Klugkist from the
# Dutch Organization for Scientific research (NWO 452-12-010).
# ----------------------------------------------------------

library(Rcpp)
library(BH)

source("DataAnalysis/Gibbs/DW.R")
source("DataAnalysis/MH/VMMH.R")
source("DataAnalysis/Rejection/FM.R")
source("DataAnalysis/describeCirc.R")
source("Simulation/vonMises.R")
Rcpp::sourceCpp("DataAnalysis/VenterMode.cpp")
Rcpp::sourceCpp("Data/rvmc.cpp")


# Perform simulations.
simulateVM <- function (nsim, n, kappa, J, meandif, Q,
                        burn, lag=1, FUN, prob = .95, printsim=FALSE,
                        wd = getwd(), ...) {
  # FUNCTION simulateVM ----------------------------------------------------
  # nsim: The number of datasets to read in. These must have been generated
  #       beforehand.
  # n, kappa, J, meandif: Properties of the datasets to be read in for analysis.
  #                       Respectively the sample size, concentration, number
  #                       of groups and difference between the group means.
  # Q: The desired number of iterations to run the chosen MCMC method.
  # burn: Number of iterations to discard as burn-in.
  # lag: Number representing a thinning factor.
  #      Only 1/lag iterations will be saved.
  # prob: We compute a prob*100% credible interval within. For example, with
  #       prob=.95, the symmetric interval used is (.025, .975).
  # printsim: Whether to print the number of the current iteration.
  # wd: Current working directory.
  # Returns: A list, with a matrix of mean directions, a vector of
  #          concentrations, and a list, 'spec' of additional values.
  # ------------------------------------------------------------------------


  # The probabilities used later for CCI's.
  probs <- c((1-prob)/2, 1-(1-prob)/2)

  # Create room for output.
  # Results to be saved for each group
  sample_mean_directions <- mu_post_mean <- sample_approx_kaps <-
    LB_mean_CCI <- UB_mean_CCI <- true_mean_in_CCI <- matrix(NA, nrow=nsim, ncol=J)

  colnames(sample_mean_directions) <- paste0("Sample mean gr ", 1:J)
  colnames(mu_post_mean)           <- paste0("Posterior mean gr ", 1:J)
  colnames(sample_approx_kaps)     <- paste0("Sample approx kap gr ", 1:J)
  colnames(LB_mean_CCI)            <- paste0("LB ", prob*100, "% Mean CCI gr ", 1:J)
  colnames(UB_mean_CCI)            <- paste0("UB ", prob*100, "% Mean CCI gr ", 1:J)
  colnames(true_mean_in_CCI)       <- paste0("True mean in ", prob*100, "% CCI gr ", 1:J)

  # Results to be saved in a nsim*2 matrix.
  kap_post_HDI <- kap_post_CI <- matrix(nrow=nsim, ncol=2)
  colnames(kap_post_CI)  <- paste0("Kap post ", prob*100, "% CI ", c("LB", "UB"))
  colnames(kap_post_HDI) <- paste0("Kap post ", prob*100, "% HDI ", c("LB", "UB"))

  # Results to be saved in a vector.
  ct <- crashed <- centred_approxKappaML <- true_kap_in_CI <- true_kap_in_HDI  <-
    R_n_post <- kap_post_mean <- kap_post_median <-  kap_post_mode.1 <- kap_post_mode.05 <-
    kap_post_mode.01 <- acceptance <- rep(NA, nsim)

  # The true means.
  true_means <- (meandif*1:J)%%(2*pi)

  # Under the uninformative prior, all prior properties are zero.
  n0 <- rep(0, J)

  for (i in 1:nsim){
    if (printsim) cat(i, "-", sep="")

    res <- list(mu =  matrix(NA, nrow=Q, ncol=J), kappa = rep(NA, Q), spec=list())

    # Read in data data.
    readfilename <- paste0(wd, "/Data/Datasets/Datasets_",
                           "J=", J, "_n=", n, "_kap=", kappa, "/nr", i, ".csv")
    th <- read.table(readfilename, sep=",")

    # Data properties.
    sample_mean_directions[i, ] <- sapply(th, meanDir)
    sample_approx_kaps[i, ]     <- sapply(th, approxKappaML)

    # Combine centered data and estimate kappa.
    centth <- sapply(1:J, function(j) th[[j]]-sample_mean_directions[i, j])
    centred_approxKappaML[i] <- approxKappaML(as.numeric(centth))


    # Run the MCMC function. By default, we take uninformative priors.
    # try() will catch the errors, system.time() will save duration.
    try({
      ct[i] <- system.time(res <- FUN(th=th, R_0=n0, mu_0=n0, c=n0, Q=Q,
                                      burn=burn, lag=lag, ...) )[3]
    },silent=TRUE)


    # If the run has crashed, we'll go on to the next attempt. It may have
    # crashed by producing NA's somewhere along the chain, for example due to
    # over/underflow, in which case the mean of the iteration means will be NA
    # as well. Another option could be that an error is produced by one of the
    # functions, in which case the try statement will catch the error and ct[i]
    # will remain NA. If either of these happens, we will not compute any
    # other parameters.
    crashed[i] <- is.na(ct[i]) || any(is.na(res$mu))

    if(!crashed[i]) {

      if(!is.null(res$spec$acceptance)) {
        acceptance[i] <- res$spec$acceptance
      } else {
        acceptance[i] <- 1
      }

      R_n_post[i] <- res$spec$R_t

      ### MU
      # Mean direction of the posterior mu's.
      mu_post_mean[i, ] <- apply(res$mu, 2, meanDir)

      # Because of the circular nature of the sample space of mean directions,
      # care must be taken here. When LB < UB, the true mean is in the CI if mu
      # > LB [AND] mu < UB. However, if LB > UB (which happens when zero degrees
      # is in the CI), the true mean is in the CI when mu > LB [OR] mu < UB.
      for (j in 1:J) {

        # The circularQuantile function from describeCirc.R ensures that the
        # correct circular quantile is obtained, which differs from the linear quantile.
        mean_CCI <- circularQuantile(res$mu[,j], probs=probs, na.rm=TRUE)

        LBij <- mean_CCI[1] %% (2*pi)
        UBij <- mean_CCI[2] %% (2*pi)
        if (LBij < UBij) {
          true_mean_in_CCI[i,j] <-  true_means[j] > LBij & true_means[j] < UBij
        } else {
          true_mean_in_CCI[i,j] <-  true_means[j] > LBij | true_means[j] < UBij
        }
        LB_mean_CCI[i,j] <- LBij
        UB_mean_CCI[i,j] <- UBij
      }

      ### KAPPA
      # Three different ways to estimate the posterior value for kappa.
      kap_post_mean[i]        <- mean(res$kappa)
      kap_post_median[i]      <- median(res$kappa)

      # The mode for continuous data requires a tuning parameter, so we try two
      # values. hmode implements the "Venter" mode in Rcpp.
      kap_post_mode.1[i]   <- hmode(res$kappa, cip=.1)
      kap_post_mode.05[i]  <- hmode(res$kappa, cip=.05)
      kap_post_mode.01[i]  <- hmode(res$kappa, cip=.01)

      # Normal confidence interval (based on the median/quantile).
      kap_post_CI[i, ]     <- quantile(res$kappa, probs = probs)

      # The Highest Density Interval (HDI), a CI for the mode.
      kap_post_HDI[i, ]    <- hmodeci(res$kappa, cip=prob)
      true_kap_in_CI[i]    <- (kappa > kap_post_CI[i, 1] &
                                 kappa < kap_post_CI[i, 2])
      true_kap_in_HDI[i]   <- (kappa > kap_post_HDI[i, 1] &
                                 kappa < kap_post_HDI[i, 2])



    }
  }

  result <- cbind("Simnr" = 1:nsim,
                  matrix(true_means, ncol=J, nrow=nsim, byrow=TRUE,
                         dimnames=list(NULL,paste0("True mean gr ", 1:J))),
                  sample_mean_directions,
                  mu_post_mean,
                  LB_mean_CCI,
                  UB_mean_CCI,
                  true_mean_in_CCI,
                  "True Kappa"=rep(kappa, nsim),
                  sample_approx_kaps,
                  "Approx Cent. KappaML" = centred_approxKappaML,
                  "Kap post mean" = kap_post_mean,
                  "Kap post median" = kap_post_median,
                  "Kap post mode.1" = kap_post_mode.1,
                  "Kap post mode.05" = kap_post_mode.05,
                  "Kap post mode.01" = kap_post_mode.01,
                  kap_post_CI,
                  kap_post_HDI,
                  "Kap true in CI" = true_kap_in_CI,
                  "Kap true in HDI" =true_kap_in_HDI,
                  "Post total R_n" = R_n_post,
                  "Acc. ratio" = acceptance,
                  "Crashed" = crashed,
                  "Comp. time" = ct)
  result
}



