#include <Rcpp.h>
#include <iostream>
#include <Data/rvmc.cpp>
#include <boost/math/special_functions/bessel.hpp>
using namespace Rcpp;

  // Obtain the log-likelihood for the posterior of the von Mises distribution.
double dVonMisesLogPosterior (double mu, double kp, double mu_n, double R_n, int m) {
    return kp * R_n * cos(mu - mu_n) - m * log(boost::math::cyl_bessel_i(0, kp));
}


// The logarithm of the chi-squared pdf, which is used as the proposal.
double dlogchisqc (double x, double df) {
  double halfdf = df/2;
  return (halfdf - 1) * log(x) - x / 2 - halfdf * log(2) - lgamma(halfdf);
}

()

//  MH sampler
// [[Rcpp::export]]
  Rcpp::List VMMHC(double kp_start,
                   NumericVector mu_n, NumericVector R_n, NumericVector m,
                   double R_t, int Qb, int lag) {

  int J = mu_n.size();

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

    // Calculate the log of the Metropolis-Hastings Ratio.
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

  NumericMatrix out(Qb, 1 + J);
  out( _, 0) = kp;
  for (int mui = 0; mui < J; mui++) {
    out( _, mui + 1) = mu( _, mui);
  }


  return Rcpp::List::create(Rcpp::Named("sam") = out,
                            Rcpp::Named("acc") = Rcpp::wrap(accepted));
}

