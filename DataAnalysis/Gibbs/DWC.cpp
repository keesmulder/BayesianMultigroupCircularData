/*
----------------------------------------------------------
DWC.cpp
MCMC-Sampling algorithm for circular data, following Damien & Walker (1999).

Kees Tim Mulder
Last updated: November 2014

This work was supported by a Vidi grant awarded to I. Klugkist from the
Dutch Organization for Scientific research (NWO 452-12-010).
----------------------------------------------------------
*/

#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>


// [[Rcpp::export]]
NumericMatrix DWC(double start_mu, double start_w, double start_kappa,
                  NumericVector lambda,
                  NumericVector mu_n, NumericVector R, double R_n,
                  int m, int Z, int Q, int lag) {
  /* FUNCTION DWC -------------------------------------------
  Generates samples from the posterior of k von Mises distributions, each with
  a separate mean, but with one common concentration kappa.
  ------------------------------------------------------------ */

  // Save the number of groups.
  int k = mu_n.size();

  // Create empty output matrices.
  NumericMatrix mu(Q, k);
  NumericVector w(Q);
  NumericVector kappa(Q);

  // Set tmu, the variable that saves the current value of mu, which is only
  // copied if the current iteration is not thinned out.
  NumericVector tmu(k);
  tmu = start_mu;

  // Other variable initializations.
  double logtau, Rsum, g, mu_min, mu_max, M, U1, v_n, Rjdis,
         kappa_min, F_k, N_k, N, acosg, odmmo, tkappa;
  long double tw, w_min, U2;
  tw     = start_w;
  tkappa = start_kappa;

  // odmmo is always one divided by m-1.
  odmmo = (1.0/(m-1));

  // Compute number of iterations, taking lag into account.
  int Qbylag = Q * lag;
  int idvlag = 1;

  // Algorithm as described in the paper.
  for (int i = 0; i < Qbylag; i++){

    // Variables whose values are used multiple times.
    logtau = log (runif(1)[0]);

    //// Sample a value for the mean for each group. ////
    for (int j = 0; j < k; j++) {
      Rjdis = R[j] * (1 + cos(tmu[j] - mu_n[j]));
      g     = (logtau/(R[j] * tkappa)) + (Rjdis/R[j]) - 1;
      if (g < -1) {
        g = -1;
      }
      acosg = acos (g);

      mu_min = mu_n[j] - acosg;
      mu_max = mu_n[j] + acosg;
      tmu[j] = runif (1, mu_min, mu_max)[0];
    }

    Rsum = sum(R * (1 + cos(tmu - mu_n)));

    //// Sample a new value for w. ////
    M = tw + rexp (1, boost::math::cyl_bessel_i(0, tkappa) - 1 )[0];
    w_min = pow (runif(1,0,1)[0], odmmo) * tw;

    U1 = runif (1, 0, 1 - exp(w_min - M))[0];
    tw = -log(1-U1) + w_min;


    //// Sample a new value for kappa. ////
    v_n = (logtau/Rsum) + tkappa;

    kappa_min = 0;
    if (v_n > 0) {
      kappa_min = v_n;
    }

    F_k = rexp(1, tw*lambda[0]*pow(tkappa, 2))[0];
    N   = tkappa * pow(1 + F_k, 0.5);

    // Choosing the correct N_k.
    for (int k_att = 2; k_att <= Z; k_att++) {
      F_k = rexp(1, tw*lambda[k_att-1]*pow(tkappa, 2*k_att))[0];
      N_k = tkappa * pow(1 + F_k, 1/(2.0*k_att));

      if ( (N_k < N)) {
        N = N_k;
      }
    }

    // The new value for kappa.
    U2     = runif(1, 0, 1 - exp(- R_n * (N - kappa_min)))[0] ;
    tkappa = -(log(1-U2)/R_n) + kappa_min;


    // For non-thinned out iterations, save the current values.
    if (i % lag == 0) {
     idvlag = i/lag;
     mu(idvlag, _) = tmu;
     w[idvlag]     = tw;
     kappa[idvlag] = tkappa;
    }

  }

  // Gather output.
  NumericMatrix out(Q, 2 + k);
  out( _, 0) = w;
  out( _, 1) = kappa;
  for (int mui = 0; mui < k; mui++) {
    out( _, mui + 2) = mu( _, mui);
  }

  return out;
}


