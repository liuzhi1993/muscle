#include <Rcpp.h>

using namespace Rcpp;

namespace muscle_pst {
List MUSCLE(NumericVector&, NumericVector&, double, bool, bool);
List DMUSCLE(NumericVector&, NumericVector&, double, int, bool, bool);
List MMUSCLE(NumericVector&, NumericMatrix&, NumericVector&, bool, bool);
List MUSCLE_FULL(NumericVector&, NumericVector&, double, bool, bool);
List DMUSCLE_FULL(NumericVector&, NumericVector&, double, int, bool, bool);
}

namespace muscle_prt {
List MUSCLE(NumericVector&, NumericVector&, double, bool, bool);
List DMUSCLE(NumericVector&, NumericVector&, double, int, bool, bool);
List MMUSCLE(NumericVector&, NumericMatrix&, NumericVector&, bool, bool);
List MUSCLE_FULL(NumericVector&, NumericVector&, double, bool, bool);
List DMUSCLE_FULL(NumericVector&, NumericVector&, double, int, bool, bool);
}

namespace muscle_original {
List MUSCLE(NumericVector&, NumericVector&, double, bool, bool);
List DMUSCLE(NumericVector&, NumericVector&, double, int, bool, bool);
List MMUSCLE(NumericVector&, NumericMatrix&, NumericVector&, bool, bool);
List MUSCLE_FULL(NumericVector&, NumericVector&, double, bool, bool);
List DMUSCLE_FULL(NumericVector&, NumericVector&, double, int, bool, bool);
}

namespace {
void check_backend(int backend) {
  if(backend < 0 || backend > 2) {
    stop("Unknown MUSCLE backend. Use 0 (PRT), 1 (PST), or 2 (original).");
  }
}
}

List muscle_dispatch(NumericVector& Y, NumericVector& q, double beta,
                     bool test, bool details, int backend) {
  check_backend(backend);
  if(backend == 0) return muscle_prt::MUSCLE(Y, q, beta, test, details);
  if(backend == 1) return muscle_pst::MUSCLE(Y, q, beta, test, details);
  return muscle_original::MUSCLE(Y, q, beta, test, details);
}

List dmuscle_dispatch(NumericVector& Y, NumericVector& q, double beta, int lag,
                      bool test, bool details, int backend) {
  check_backend(backend);
  if(backend == 0) return muscle_prt::DMUSCLE(Y, q, beta, lag, test, details);
  if(backend == 1) return muscle_pst::DMUSCLE(Y, q, beta, lag, test, details);
  return muscle_original::DMUSCLE(Y, q, beta, lag, test, details);
}

List mmuscle_dispatch(NumericVector& Y, NumericMatrix& q_matrix,
                      NumericVector& beta_vec, bool test, bool details,
                      int backend) {
  check_backend(backend);
  if(backend == 0) return muscle_prt::MMUSCLE(Y, q_matrix, beta_vec, test, details);
  if(backend == 1) return muscle_pst::MMUSCLE(Y, q_matrix, beta_vec, test, details);
  return muscle_original::MMUSCLE(Y, q_matrix, beta_vec, test, details);
}

List muscle_full_dispatch(NumericVector& Y, NumericVector& q, double beta,
                          bool test, bool details, int backend) {
  check_backend(backend);
  if(backend == 0) return muscle_prt::MUSCLE_FULL(Y, q, beta, test, details);
  if(backend == 1) return muscle_pst::MUSCLE_FULL(Y, q, beta, test, details);
  return muscle_original::MUSCLE_FULL(Y, q, beta, test, details);
}

List dmuscle_full_dispatch(NumericVector& Y, NumericVector& q, double beta,
                           int lag, bool test, bool details, int backend) {
  check_backend(backend);
  if(backend == 0) return muscle_prt::DMUSCLE_FULL(Y, q, beta, lag, test, details);
  if(backend == 1) return muscle_pst::DMUSCLE_FULL(Y, q, beta, lag, test, details);
  return muscle_original::DMUSCLE_FULL(Y, q, beta, lag, test, details);
}
