#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <cmath>

// Vectorized Poisson PMF (like R's dpois)
// [[Rcpp::export]]
std::vector<double> mydpois(const std::vector<double>& x,
                            const std::vector<double>& lambda,
                            bool log = false) {
  if (x.empty() || lambda.empty()) {
    throw std::invalid_argument("x and lambda must not be empty.");
  }

  // Recycle vectors to the maximum length, R-style
  size_t n = std::max(x.size(), lambda.size());
  std::vector<double> result(n);

  for (size_t i = 0; i < n; ++i) {
    double xi = x[i % x.size()];
    double lam = lambda[i % lambda.size()];

    double intpart;
    bool not_integer = std::modf(xi, &intpart) != 0.0;

    if (lam <= 0.0 || xi < 0 || not_integer) {
      result[i] = log ? -INFINITY : 0.0;
      if (not_integer) {
        std::cerr << "Warning: non-integer x = " << xi << std::endl;
      }
    } else {
      // Safe integer conversion
      int k = static_cast<int>(xi);
      double logp = -lam + k * std::log(lam) - std::lgamma(k + 1.0);
      result[i] = log ? logp : std::exp(logp);
    }
  }

  return result;
}

