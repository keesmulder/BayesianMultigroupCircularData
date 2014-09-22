#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;


// Calculate the mode of a sample of data as in Venter (1967), to obtain the
// HPD, for example.

// Find the highest density point.

// [[Rcpp::export]]
double hmode(NumericVector x, double cip) {

  int n, cil, chiv;
  double ln, M;

  n = x.size();
  NumericVector sx = clone(x);
  std::sort(sx.begin(), sx.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Size of the currently smallest interval.
  ln = sx[cil]-sx[0];

  for (int i=0; i < (n-cil); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (sx[i+cil]-sx[i])) {
      ln = (sx[i+cil]-sx[i]);
      chiv = i;
    }
  }

  M = (sx[chiv+cil]+sx[chiv])/2;

  return M;
}



// Find the highest density interval.

// [[Rcpp::export]]
NumericVector hmodeci(NumericVector x, double cip) {

  int n, cil, chiv;
  double ln;

  n = x.size();
  NumericVector sx = clone(x);
  std::sort(sx.begin(), sx.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Size of the currently smallest interval.
  ln = sx[cil]-sx[0];

  for (int i=0; i < (n-cil); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (sx[i+cil]-sx[i])) {
      ln = (sx[i+cil]-sx[i]);
      chiv = i;
    }
  }

  NumericVector M(2);
  M[0] = sx[chiv];
  M[1] = sx[chiv+cil];

  return M;
}
