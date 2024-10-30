// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// kernal
double kernal(double dx);
RcppExport SEXP _SSIndex3_kernal(SEXP dxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type dx(dxSEXP);
    rcpp_result_gen = Rcpp::wrap(kernal(dx));
    return rcpp_result_gen;
END_RCPP
}
// M2
double M2(const arma::vec& beta, const arma::mat& Z, const arma::vec& T2, const arma::vec& C2, const arma::vec& Y2, double h, double ht);
RcppExport SEXP _SSIndex3_M2(SEXP betaSEXP, SEXP ZSEXP, SEXP T2SEXP, SEXP C2SEXP, SEXP Y2SEXP, SEXP hSEXP, SEXP htSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type T2(T2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type C2(C2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y2(Y2SEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type ht(htSEXP);
    rcpp_result_gen = Rcpp::wrap(M2(beta, Z, T2, C2, Y2, h, ht));
    return rcpp_result_gen;
END_RCPP
}
// shapeFun3
double shapeFun3(int n, arma::vec& m, arma::vec& midx, arma::vec& tij, arma::vec& yi, arma::mat& xb, double x, double t, double h, arma::vec& w, arma::vec& med, double result);
RcppExport SEXP _SSIndex3_shapeFun3(SEXP nSEXP, SEXP mSEXP, SEXP midxSEXP, SEXP tijSEXP, SEXP yiSEXP, SEXP xbSEXP, SEXP xSEXP, SEXP tSEXP, SEXP hSEXP, SEXP wSEXP, SEXP medSEXP, SEXP resultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type midx(midxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tij(tijSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xb(xbSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type med(medSEXP);
    Rcpp::traits::input_parameter< double >::type result(resultSEXP);
    rcpp_result_gen = Rcpp::wrap(shapeFun3(n, m, midx, tij, yi, xb, x, t, h, w, med, result));
    return rcpp_result_gen;
END_RCPP
}
// size
double size(int n, arma::mat& xr, arma::vec& mFhat, arma::vec& w);
RcppExport SEXP _SSIndex3_size(SEXP nSEXP, SEXP xrSEXP, SEXP mFhatSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mFhat(mFhatSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(size(n, xr, mFhat, w));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SSIndex3_kernal", (DL_FUNC) &_SSIndex3_kernal, 1},
    {"_SSIndex3_M2", (DL_FUNC) &_SSIndex3_M2, 7},
    {"_SSIndex3_shapeFun3", (DL_FUNC) &_SSIndex3_shapeFun3, 12},
    {"_SSIndex3_size", (DL_FUNC) &_SSIndex3_size, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SSIndex3(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
