// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// x1
double x1(double q, double b);
RcppExport SEXP _MUSCLE_x1(SEXP qSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(x1(q, b));
    return rcpp_result_gen;
END_RCPP
}
// x2
double x2(double q, double b);
RcppExport SEXP _MUSCLE_x2(SEXP qSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(x2(q, b));
    return rcpp_result_gen;
END_RCPP
}
// MUSCLE
List MUSCLE(NumericVector& Y, NumericVector& q, double b, bool test);
RcppExport SEXP _MUSCLE_MUSCLE(SEXP YSEXP, SEXP qSEXP, SEXP bSEXP, SEXP testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type test(testSEXP);
    rcpp_result_gen = Rcpp::wrap(MUSCLE(Y, q, b, test));
    return rcpp_result_gen;
END_RCPP
}
// DMUSCLE
List DMUSCLE(NumericVector& Y, NumericVector& q, double b, int lag, bool test);
RcppExport SEXP _MUSCLE_DMUSCLE(SEXP YSEXP, SEXP qSEXP, SEXP bSEXP, SEXP lagSEXP, SEXP testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< bool >::type test(testSEXP);
    rcpp_result_gen = Rcpp::wrap(DMUSCLE(Y, q, b, lag, test));
    return rcpp_result_gen;
END_RCPP
}
// MMUSCLE
List MMUSCLE(NumericVector& Y, NumericMatrix& q, NumericVector& b, bool test);
RcppExport SEXP _MUSCLE_MMUSCLE(SEXP YSEXP, SEXP qSEXP, SEXP bSEXP, SEXP testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type test(testSEXP);
    rcpp_result_gen = Rcpp::wrap(MMUSCLE(Y, q, b, test));
    return rcpp_result_gen;
END_RCPP
}
// MUSCLE_FULL
List MUSCLE_FULL(NumericVector& Y, NumericVector& q, double b, bool test);
RcppExport SEXP _MUSCLE_MUSCLE_FULL(SEXP YSEXP, SEXP qSEXP, SEXP bSEXP, SEXP testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type test(testSEXP);
    rcpp_result_gen = Rcpp::wrap(MUSCLE_FULL(Y, q, b, test));
    return rcpp_result_gen;
END_RCPP
}
// DMUSCLE_FULL
List DMUSCLE_FULL(NumericVector& Y, NumericVector& q, double b, int lag, bool test);
RcppExport SEXP _MUSCLE_DMUSCLE_FULL(SEXP YSEXP, SEXP qSEXP, SEXP bSEXP, SEXP lagSEXP, SEXP testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< bool >::type test(testSEXP);
    rcpp_result_gen = Rcpp::wrap(DMUSCLE_FULL(Y, q, b, lag, test));
    return rcpp_result_gen;
END_RCPP
}
// simulQuantile_MUSCLE
NumericVector simulQuantile_MUSCLE(double p, int n, double b);
RcppExport SEXP _MUSCLE_simulQuantile_MUSCLE(SEXP pSEXP, SEXP nSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(simulQuantile_MUSCLE(p, n, b));
    return rcpp_result_gen;
END_RCPP
}
// logg
NumericVector logg(NumericVector& x);
RcppExport SEXP _MUSCLE_logg(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logg(x));
    return rcpp_result_gen;
END_RCPP
}
// simulQuantile_DMUSCLE
NumericVector simulQuantile_DMUSCLE(NumericVector& X, NumericVector& ACF, int n, double b);
RcppExport SEXP _MUSCLE_simulQuantile_DMUSCLE(SEXP XSEXP, SEXP ACFSEXP, SEXP nSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type ACF(ACFSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(simulQuantile_DMUSCLE(X, ACF, n, b));
    return rcpp_result_gen;
END_RCPP
}
// simulQuantile_DMUSCLE2
NumericVector simulQuantile_DMUSCLE2(double p, int n, NumericVector& ACF, NumericMatrix& L);
RcppExport SEXP _MUSCLE_simulQuantile_DMUSCLE2(SEXP pSEXP, SEXP nSEXP, SEXP ACFSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type ACF(ACFSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(simulQuantile_DMUSCLE2(p, n, ACF, L));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MUSCLE_x1", (DL_FUNC) &_MUSCLE_x1, 2},
    {"_MUSCLE_x2", (DL_FUNC) &_MUSCLE_x2, 2},
    {"_MUSCLE_MUSCLE", (DL_FUNC) &_MUSCLE_MUSCLE, 4},
    {"_MUSCLE_DMUSCLE", (DL_FUNC) &_MUSCLE_DMUSCLE, 5},
    {"_MUSCLE_MMUSCLE", (DL_FUNC) &_MUSCLE_MMUSCLE, 4},
    {"_MUSCLE_MUSCLE_FULL", (DL_FUNC) &_MUSCLE_MUSCLE_FULL, 4},
    {"_MUSCLE_DMUSCLE_FULL", (DL_FUNC) &_MUSCLE_DMUSCLE_FULL, 5},
    {"_MUSCLE_simulQuantile_MUSCLE", (DL_FUNC) &_MUSCLE_simulQuantile_MUSCLE, 3},
    {"_MUSCLE_logg", (DL_FUNC) &_MUSCLE_logg, 1},
    {"_MUSCLE_simulQuantile_DMUSCLE", (DL_FUNC) &_MUSCLE_simulQuantile_DMUSCLE, 4},
    {"_MUSCLE_simulQuantile_DMUSCLE2", (DL_FUNC) &_MUSCLE_simulQuantile_DMUSCLE2, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_MUSCLE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
