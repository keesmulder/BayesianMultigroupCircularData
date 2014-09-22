#include <Rcpp.h>
#include <iostream>
#include <boost/math/special_functions/bessel.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
  NumericMatrix DWC(double start_mu, double start_w, double start_kappa,
                  NumericVector lambda,
                  NumericVector mu_n, NumericVector R, double R_n,
                  int m, int N_k_attempts, int Q, int lag) {

  int k = mu_n.size();

  NumericMatrix mu(Q, k);
  NumericVector w(Q);
  NumericVector kappa(Q);

  // The variable that saves the current value of mu, which
  // is only copied if the current iteration is not thinned
  // out.
  NumericVector tmu(k);


  //for (int mui = 0; mui < k; mui++) {
  //  mu(0, mui) = start_mu;
  //}
  tmu = start_mu;

  // Variable definitions.
  double logtau, Rsum, g, mu_min, mu_max, M, U1, v_n, kappa_min, F_k, N_k, N, acosg, odmmo, tkappa;
  long double tw, w_min, U2;
  tw     = start_w;
  tkappa = start_kappa;

  // One divided by m-1.
  odmmo = (1.0/(m-1));

  // Number of iterations given lag.
  int Qbylag = Q * lag;
  int idvlag = 1;

  for (int i = 0; i < Qbylag; i++){

    // Algorithm as described in the paper.

    logtau = log (runif(1)[0]);
    Rsum = sum(R * (1 + cos(tmu - mu_n)));

    // Sample means
    g = (logtau/(R_n * tkappa)) + (Rsum/R_n) - 1;
    if (g < -1) {
      g = -1;
    }
    acosg = acos (g);

    for (int j = 0; j < k; j++) {
      mu_min = mu_n[j] - acosg;
      mu_max = mu_n[j] + acosg;
      tmu[j] = runif (1, mu_min, mu_max)[0];
    }

    M = tw + rexp (1, boost::math::cyl_bessel_i(0, tkappa) - 1 )[0];
    w_min = pow (runif(1,0,1)[0], odmmo) * tw;

    U1    = runif (1, 0, 1 - exp(w_min - M))[0];
    tw = -log(1-U1) + w_min;

    v_n = (logtau/Rsum) + tkappa;

    kappa_min = 0;
    if (v_n > 0) {
      kappa_min = v_n;
    }

    F_k = rexp(1, tw*lambda[0]*pow(tkappa, 2))[0];
    N   = tkappa * pow(1 + F_k, 0.5);

    // Choosing the correct N_k.
    for (int k_att = 2; k_att <= N_k_attempts; k_att++) {
      F_k = rexp(1, tw*lambda[k_att-1]*pow(tkappa, 2*k_att))[0];
      N_k = tkappa * pow(1 + F_k, 1/(2.0*k_att));

      if ( (N_k < N)) {
        N = N_k;
      }
    }

    U2        = runif(1, 0, 1 - exp(- R_n * (N - kappa_min)))[0] ;

    // The new value for kappa.
    tkappa = -(log(1-U2)/R_n) + kappa_min;

    //  For non-thinned out iterations, save the current value.
    if (i % lag == 0) {
     idvlag = i/lag;
     mu(idvlag, _) = tmu;
     w[idvlag] = tw;
     kappa[idvlag] = tkappa;
    }

  }

  // Output.
  NumericMatrix out(Q, 2 + k);
  out( _, 0) = w;
  out( _, 1) = kappa;
  for (int mui = 0; mui < k; mui++) {
    out( _, mui + 2) = mu( _, mui);
  }

  return out;
}


