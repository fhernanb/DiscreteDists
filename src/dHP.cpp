#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <cmath>
#include <complex>

//' Function to obtain F11 with C++.
//' @param gamma numeric value for gamma.
//' @param lambda numeric value for lambda.
//' @param maxiter_series numeric value.
//' @param tol numeric value.
//' @export
//' @return returns the F11 value.
// [[Rcpp::export]]
 double f11_cpp(double gamma, double lambda, int maxiter_series = 10000, double tol = 1.0e-10) {
   double fac  = 1.0;
   double temp = 1.0;
   double L    = gamma;
   double series = temp;
   double f11;

   for (int n = 1; n <= maxiter_series; ++n) {
     fac = fac * lambda / L;
     series = temp + fac;

     if (std::abs(series - temp) < tol) {
       f11 = series;  // Assuming series is already real in this context
       return f11;
     }

     temp = series;
     L += 1;
   }

   f11 = series;  // Assuming series is already real in this context
   return f11;
 }

//' Function to obtain the dHYPERPO for a single value x
//' @param x numeric value for x.
//' @param mu numeric value for nu.
//' @param sigma numeric value for sigma.
//' @param log logical value for log.
//' @export
//' @return returns the pmf for a single value x.
// [[Rcpp::export]]
 double dHYPERPO_single(double x, double mu=1, double sigma=1, bool log=false) {
   if (sigma <= 0 || mu <= 0) {
     throw std::runtime_error("parameter sigma and mu must be positive!");
   }

   double res;

   if (x < 0) {
     res = std::log(0);
   }
   else {
     double p1 = x * std::log(mu) - lgamma(sigma + x) + lgamma(sigma);
     double temp_f11 = f11_cpp(sigma, mu);
     double p2 = std::log(temp_f11);
     res = p1 - p2;
   }

   if (log) {
     return res;
   } else {
     return std::exp(res);
   }
 }

//' Function to obtain the dHYPERPO for a vector x
//' @param x numeric value for x.
//' @param mu numeric value for nu.
//' @param sigma numeric value for sigma.
//' @param log logical value for log.
//' @export
//' @return returns the pmf for a vector.
// [[Rcpp::export]]
 NumericVector dHYPERPO_vec(NumericVector x, NumericVector mu,
                            NumericVector sigma, LogicalVector log) {
   int n = x.size();
   NumericVector out(n);

   for(int i = 0; i < n; ++i) {
     out[i] = dHYPERPO_single(x[i], mu[i], sigma[i], log[i]);
   }
   return out;
 }

