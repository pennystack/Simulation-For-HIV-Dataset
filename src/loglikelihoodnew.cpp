#include <iostream>
#include <Rcpp.h>
#include <math.h> 
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
double cpp_p(int i, int j, NumericVector aij, NumericVector bij, double t, int nstates) {
  double out = t * aij[(j - 1) * nstates + i - 1] + bij[(j - 1) * nstates + i - 1];
  //if (out < 1e-100) return(1e-100);
  //if (out > 1e5) return(1);
  
  return std::exp(out) / (1 + std::exp(out));
}

// [[Rcpp::export]]
double cpp_f(int i, int j, NumericVector vij, NumericVector sij, double x, int nstates) {
  double out = (vij[(j - 1) * nstates + i - 1] / sij[(j - 1) * nstates + i - 1]) * std::pow(x / sij[(j - 1) * nstates + i - 1], vij[(j - 1) * nstates + i - 1] - 1) * std::exp(-std::pow(x / sij[(j - 1) * nstates + i - 1], vij[(j - 1) * nstates + i - 1]));
  if (out < 1e-100) out = 1e-100;

  return out;
}

// [[Rcpp::export]]
double cpp_F(int i, int j, NumericVector vij, NumericVector sij, double x, int nstates) {
  double out = 1 - std::exp(-std::pow(x / sij[(j - 1) * nstates + i - 1], vij[(j - 1) * nstates + i - 1]));
  if (out < 1e-100) out = 1e-100;
  if (out > 1) out = 1;
  
  return out;
}

// [[Rcpp::export]]
double cpp_H(int i, double t, double x, NumericVector aij, NumericVector bij, NumericVector vij, NumericVector sij, int nstates) {
  
  double out = 0.0;
  for (int j = 1; j <= nstates; j++) {
    if (j != i) {
      out += cpp_p(i, j, aij, bij, t, nstates) * cpp_F(i, j, vij, sij, x, nstates);
    }
  }
  
  if (out < 1e-100) out = 1e-100;
  if (out > 1) out = 1;
  
  return out;
}

// [[Rcpp::export]]
double cpp_S(int i, double t, double x, NumericVector aij, NumericVector bij, NumericVector vij, NumericVector sij, int nstates) {
  
  double out = 1 - cpp_H(i, t, x, aij, bij, vij, sij, nstates);
  
  if (out < 1e-100) out = 1e-100;
  if (out > 1) out = 1;
  
  return out;
}

// [[Rcpp::export]]
double cpp_loglikelihood(DataFrame obstimes, NumericVector aij, NumericVector bij, NumericVector vij, NumericVector sij, int nstates) {
  NumericVector patients = obstimes["PATIENT"] ;
  NumericVector state = obstimes["state"] ;
  NumericVector obstime = obstimes["obstime"] ;
  NumericVector deltaobstime = obstimes["deltaobstime"] ;
  
  double out = 0.;
  
  int n = patients.length() ;
  
  for (int i = 0; i < n; i++) {
    if (patients[i] == patients[i + 1] and state[i] != state[i + 1]) {
      out += (log(cpp_p(state[i], state[i + 1], aij, bij, obstime[i], nstates)) + log(cpp_f(state[i], state[i + 1], vij, sij, deltaobstime[i], nstates))) ;
    }
    if (patients[i] == patients[i + 1] and state[i] == state[i + 1]) {
      out += log(cpp_S(state[i], obstime[i], deltaobstime[i], aij, bij, vij, sij, nstates)) ;
    }
  }
  
  return out;
}


// [[Rcpp::export]]
double cpp_c(int i, int j, NumericVector aij, NumericVector bij, NumericVector vij, NumericVector sij, double t, double x, int nstates) {
  double out = cpp_p(i, j, aij, bij, t, nstates) * cpp_f(i, j, vij, sij, x, nstates) ;
  if (out < 1e-100) out = 1e-100;
 

  return out;
}
