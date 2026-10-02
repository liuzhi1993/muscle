#include <Rcpp.h>

using namespace Rcpp;

namespace muscle_pst {
List pst_window_query(NumericVector&, int, int, int, int, int, int, int);
NumericVector simulQuantile_MUSCLE(double, int, double);
NumericVector logg(NumericVector&);
NumericVector simulQuantile_DMUSCLE(NumericVector&, NumericVector&, int);
}

List muscle_dispatch(NumericVector&, NumericVector&, double, bool, bool, int);
List dmuscle_dispatch(NumericVector&, NumericVector&, double, int, bool, bool, int);
List mmuscle_dispatch(NumericVector&, NumericMatrix&, NumericVector&, bool, bool, int);
List muscle_full_dispatch(NumericVector&, NumericVector&, double, bool, bool, int);
List dmuscle_full_dispatch(NumericVector&, NumericVector&, double, int, bool, bool, int);

RcppExport SEXP _muscle_pst_window_query(SEXP YSEXP, SEXP start_offsetSEXP, SEXP end_offsetSEXP, SEXP n_startsSEXP, SEXP orderSEXP, SEXP loSEXP, SEXP hiSEXP, SEXP paritySEXP) {
BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  NumericVector Y(YSEXP);
  return Rcpp::wrap(muscle_pst::pst_window_query(Y, as<int>(start_offsetSEXP), as<int>(end_offsetSEXP), as<int>(n_startsSEXP), as<int>(orderSEXP), as<int>(loSEXP), as<int>(hiSEXP), as<int>(paritySEXP)));
END_RCPP
}

RcppExport SEXP _muscle_MUSCLE(SEXP YSEXP, SEXP qSEXP, SEXP betaSEXP, SEXP testSEXP, SEXP detailsSEXP, SEXP backendSEXP) {
BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  NumericVector Y(YSEXP), q(qSEXP);
  return Rcpp::wrap(muscle_dispatch(Y, q, as<double>(betaSEXP), as<bool>(testSEXP), as<bool>(detailsSEXP), as<int>(backendSEXP)));
END_RCPP
}

RcppExport SEXP _muscle_DMUSCLE(SEXP YSEXP, SEXP qSEXP, SEXP betaSEXP, SEXP lagSEXP, SEXP testSEXP, SEXP detailsSEXP, SEXP backendSEXP) {
BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  NumericVector Y(YSEXP), q(qSEXP);
  return Rcpp::wrap(dmuscle_dispatch(Y, q, as<double>(betaSEXP), as<int>(lagSEXP), as<bool>(testSEXP), as<bool>(detailsSEXP), as<int>(backendSEXP)));
END_RCPP
}

RcppExport SEXP _muscle_MMUSCLE(SEXP YSEXP, SEXP q_matrixSEXP, SEXP beta_vecSEXP, SEXP testSEXP, SEXP detailsSEXP, SEXP backendSEXP) {
BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  NumericVector Y(YSEXP), beta_vec(beta_vecSEXP);
  NumericMatrix q_matrix(q_matrixSEXP);
  return Rcpp::wrap(mmuscle_dispatch(Y, q_matrix, beta_vec, as<bool>(testSEXP), as<bool>(detailsSEXP), as<int>(backendSEXP)));
END_RCPP
}

RcppExport SEXP _muscle_MUSCLE_FULL(SEXP YSEXP, SEXP qSEXP, SEXP betaSEXP, SEXP testSEXP, SEXP detailsSEXP, SEXP backendSEXP) {
BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  NumericVector Y(YSEXP), q(qSEXP);
  return Rcpp::wrap(muscle_full_dispatch(Y, q, as<double>(betaSEXP), as<bool>(testSEXP), as<bool>(detailsSEXP), as<int>(backendSEXP)));
END_RCPP
}

RcppExport SEXP _muscle_DMUSCLE_FULL(SEXP YSEXP, SEXP qSEXP, SEXP betaSEXP, SEXP lagSEXP, SEXP testSEXP, SEXP detailsSEXP, SEXP backendSEXP) {
BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  NumericVector Y(YSEXP), q(qSEXP);
  return Rcpp::wrap(dmuscle_full_dispatch(Y, q, as<double>(betaSEXP), as<int>(lagSEXP), as<bool>(testSEXP), as<bool>(detailsSEXP), as<int>(backendSEXP)));
END_RCPP
}

RcppExport SEXP _muscle_simulQuantile_MUSCLE(SEXP pSEXP, SEXP nSEXP, SEXP betaSEXP) {
BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  return Rcpp::wrap(muscle_pst::simulQuantile_MUSCLE(as<double>(pSEXP), as<int>(nSEXP), as<double>(betaSEXP)));
END_RCPP
}

RcppExport SEXP _muscle_logg(SEXP xSEXP) {
BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  NumericVector x(xSEXP);
  return Rcpp::wrap(muscle_pst::logg(x));
END_RCPP
}

RcppExport SEXP _muscle_simulQuantile_DMUSCLE(SEXP XSEXP, SEXP ACFSEXP, SEXP nSEXP) {
BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  NumericVector X(XSEXP), ACF(ACFSEXP);
  return Rcpp::wrap(muscle_pst::simulQuantile_DMUSCLE(X, ACF, as<int>(nSEXP)));
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
  {"_muscle_pst_window_query", (DL_FUNC) &_muscle_pst_window_query, 8},
  {"_muscle_MUSCLE", (DL_FUNC) &_muscle_MUSCLE, 6},
  {"_muscle_DMUSCLE", (DL_FUNC) &_muscle_DMUSCLE, 7},
  {"_muscle_MMUSCLE", (DL_FUNC) &_muscle_MMUSCLE, 6},
  {"_muscle_MUSCLE_FULL", (DL_FUNC) &_muscle_MUSCLE_FULL, 6},
  {"_muscle_DMUSCLE_FULL", (DL_FUNC) &_muscle_DMUSCLE_FULL, 7},
  {"_muscle_simulQuantile_MUSCLE", (DL_FUNC) &_muscle_simulQuantile_MUSCLE, 3},
  {"_muscle_logg", (DL_FUNC) &_muscle_logg, 1},
  {"_muscle_simulQuantile_DMUSCLE", (DL_FUNC) &_muscle_simulQuantile_DMUSCLE, 3},
  {NULL, NULL, 0}
};

extern "C" void R_init_muscle(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
