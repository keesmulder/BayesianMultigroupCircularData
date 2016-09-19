/*
----------------------------------------------------------
VMMHC.cpp
This is an RCPP implementation the MH-algorithm for multiple groups of data from
the von Mises distribution.

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

NumericVector rvmc(int n, double mu, double kp) {
  /* FUNCTION rvmc -------------------------------------------
  Generate random variates from the von Mises distribution.

  n:      The number of random variates required.
  mu:     The required mean direction, mu.
  kp:     The required concentration, kappa.

  Returns: A vector of length n containing VM random variates.
  ------------------------------------------------------------ */

  // If kappa is very small, return a circular uniform draw, as otherwise the
  // algorithm will fail.
  if (kp < .0000001) {
    return runif(n, 0, 8.0*atan(1));
  }
  
  NumericVector th(n);
  int sn;
  double a, b, r, u1, u2, u3, z, f, c;
  bool cont;

  a = 1 + sqrt(1 + 4.0 * pow(kp, 2));
  b = (a - (sqrt(2.0*a)))/(2.0*kp);
  r = (1 + pow(b,2))/(2.0*b);

  for (int i=0; i<n; i++) {

    cont = TRUE;

    do {
      u1 = runif(1, 0, 1)[0];
      u2 = runif(1, 0, 1)[0];
      u3 = runif(1, 0, 1)[0];

      // STEP 1
      z = cos(4*atan(1)*u1);
      f = (1 + r*z)/(r + z);
      c = kp * (r - f);

      // STEP 2
      if (c*(2-c) - u2 > 0) cont=FALSE;

      // STEP 3
      if (log(c/u2) + 1 - c >= 0) cont=FALSE;
    } while (cont);

    // STEP 4
    if (u3 - 0.5 > 0) {
      sn = 1;
    } else {
      sn = -1;
    }

    th[i] = fmod(sn * acos(f) + mu, 8.0*atan(1));
}

  return th;
}

// FUNCTION DEFINITIONS
// Obtain the log-likelihood for the posterior of the von Mises distribution.
// [[Rcpp::export]]
double dVonMisesLogPosterior (double mu, double kp, double mu_n, double R_n, int m) {
    return kp * R_n * cos(mu - mu_n) - m * log(boost::math::cyl_bessel_i(0, kp));
}

// Obtain the logarithm of the chi-squared pdf, which is used as the proposal.
// [[Rcpp::export]]
double dlogchisqc (double x, double df) {
  double halfdf = df/2;
  return (halfdf - 1) * log(x) - x / 2 - halfdf * log(2) - lgamma(halfdf);
}



// [[Rcpp::export]]
  Rcpp::List VMMHC(double kp_start,
                   NumericVector mu_n, NumericVector R_n, NumericVector m,
                   double R_t, int Qb, int lag) {
  /* FUNCTION VMMHC -------------------------------------------
  Generates samples from the posterior of k von Mises distributions, each with
  a separate mean, but with one common concentration kappa.
  ------------------------------------------------------------ */

  // The number of groups.
  int J = mu_n.size();

  // Create empty output matrices.
  NumericMatrix mu(Qb, J);
  NumericVector kp(Qb);

  // The variable that saves the current value of mu, which
  // is only copied if the current iteration is not thinned
  // out.
  NumericVector mu_new(J);

  double kp_cur, kp_can, post_num, post_den, prop_num, prop_den, lnMHR, u;
  int accepted = 0;

  kp_cur = kp_start;

  int Qbylag = Qb * lag;
  int idvlag = 1;

  for (int i = 0; i < Qbylag; i++){

    // Obtain a candidate for kappa
    kp_can = rchisq(1, kp_cur)[0];

    // Set values for the posterior likelihood back to 0.
    post_num = 0;
    post_den = 0;

    for (int j = 0; j < J; j++) {

      // Sample value for mu[j].
      mu_new[j] = rvmc(1, mu_n[j], R_n[j]*kp_cur)[0];

      // Calculate the numerator and denominator log-posterior of the von Mises
      // distribution, which will be used in the MH ratio.
      post_num += dVonMisesLogPosterior(mu_new[j], kp_can, mu_n[j], R_n[j], m[j]);
      post_den += dVonMisesLogPosterior(mu_new[j], kp_cur, mu_n[j], R_n[j], m[j]);
    }

    // Compute the required values of the logarithms of the proposal densities.
    prop_num = dlogchisqc(kp_cur, kp_can);
    prop_den = dlogchisqc(kp_can, kp_cur);

    // Calculate the log of the Metropolis-Hastings ratio.
    lnMHR = post_num + prop_num - prop_den - post_den;

    // Determine whether to accept the candidate.
    u = runif(1)[0];
    if (lnMHR > log(u)) {
      kp_cur = kp_can;
      accepted = accepted+1;
    }

    //    For non-thinned out iterations, save the current value.
    if (i % lag == 0) {
      idvlag = i/lag;
      mu(idvlag, _) = mu_new;
      kp[idvlag] = kp_cur;
    }
  }

  // Gather output.
  NumericMatrix out(Qb, 1 + J);
  out( _, 0) = kp;
  for (int mui = 0; mui < J; mui++) {
    out( _, mui + 1) = mu( _, mui);
  }


  return Rcpp::List::create(Rcpp::Named("sam") = out,
                            Rcpp::Named("acc") = Rcpp::wrap(accepted));
}

