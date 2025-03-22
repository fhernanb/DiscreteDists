// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// f11_cpp
double f11_cpp(double gamma, double lambda, int maxiter_series, double tol);
RcppExport SEXP _DiscreteDists_f11_cpp(SEXP gammaSEXP, SEXP lambdaSEXP, SEXP maxiter_seriesSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter_series(maxiter_seriesSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(f11_cpp(gamma, lambda, maxiter_series, tol));
    return rcpp_result_gen;
END_RCPP
}
// dHYPERPO_single
double dHYPERPO_single(double x, double mu, double sigma, bool log);
RcppExport SEXP _DiscreteDists_dHYPERPO_single(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dHYPERPO_single(x, mu, sigma, log));
    return rcpp_result_gen;
END_RCPP
}
// dHYPERPO_vec
NumericVector dHYPERPO_vec(NumericVector x, NumericVector mu, NumericVector sigma, LogicalVector log);
RcppExport SEXP _DiscreteDists_dHYPERPO_vec(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dHYPERPO_vec(x, mu, sigma, log));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DiscreteDists_f11_cpp", (DL_FUNC) &_DiscreteDists_f11_cpp, 4},
    {"_DiscreteDists_dHYPERPO_single", (DL_FUNC) &_DiscreteDists_dHYPERPO_single, 4},
    {"_DiscreteDists_dHYPERPO_vec", (DL_FUNC) &_DiscreteDists_dHYPERPO_vec, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_DiscreteDists(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
