// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// improve
arma::rowvec improve(arma::uvec& muts, double& L, arma::rowvec par);
RcppExport SEXP _segmut_improve(SEXP mutsSEXP, SEXP LSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec& >::type muts(mutsSEXP);
    Rcpp::traits::input_parameter< double& >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(improve(muts, L, par));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_segmut_improve", (DL_FUNC) &_segmut_improve, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_segmut(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}