#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <cmath>
#include <complex>


// [[Rcpp::export]]
 double dDMOLBE_single(double x, double mu=1, double sigma=1, bool log=false) {
   if (sigma <= 0 || mu <= 0) {
     throw std::runtime_error("parameter sigma and mu must be positive!");
   }

   double res;

   if (x < 0) {
     res = std::log(0);
   }
   else {
     double k1 = (1+x/mu) * std::exp(-x/mu);
     double k2 = (1+(x+1)/mu) * std::exp(-(x+1)/mu);
     double p1 = std::log(sigma) + std::log(k1 - k2);
     double p2 = -std::log(1-(1-sigma)*k1) - std::log(1-(1-sigma)*k2);
     res = p1 + p2;
   }

   if (log) {
     return res;
   } else {
     return std::exp(res);
   }
 }


// [[Rcpp::export]]
 NumericVector dDMOLBE_vec(NumericVector x, NumericVector mu,
                            NumericVector sigma, LogicalVector log) {
   int n = x.size();
   NumericVector out(n);

   for(int i = 0; i < n; ++i) {
     out[i] = dDMOLBE_single(x[i], mu[i], sigma[i], log[i]);
   }
   return out;
 }

