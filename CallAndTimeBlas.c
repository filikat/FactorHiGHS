#include "CallAndTimeBlas.h"

#include "Blas_declaration.h"
#include "DenseFact_declaration.h"
#include "timing.h"

void callAndTime_dsyrk(char uplo, char trans, int n, int k, double alpha,
                       const double* A, int lda, double beta, double* C,
                       int ldc, double* times) {
#ifdef FINEST_TIMING
  double t0 = GetTime();
#endif
  dsyrk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc);
#ifdef FINEST_TIMING
  times[kTimeDenseFact_syrk] += GetTime() - t0;
#endif
}

void callAndTime_dgemm(char transa, char transb, int m, int n, int k,
                       double alpha, const double* A, int lda, const double* B,
                       int ldb, double beta, double* C, int ldc,
                       double* times) {
#ifdef FINEST_TIMING
  double t0 = GetTime();
#endif
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C,
         &ldc);
#ifdef FINEST_TIMING
  times[kTimeDenseFact_gemm] += GetTime() - t0;
#endif
}

void callAndTime_dtrsm(char side, char uplo, char trans, char diag, int m,
                       int n, double alpha, const double* A, int lda, double* B,
                       int ldb, double* times) {
#ifdef FINEST_TIMING
  double t0 = GetTime();
#endif
  dtrsm_(&side, &uplo, &trans, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
#ifdef FINEST_TIMING
  times[kTimeDenseFact_trsm] += GetTime() - t0;
#endif
}

void callAndTime_dcopy(int n, const double* dx, int incx, double* dy, int incy,
                       double* times) {
#ifdef FINEST_TIMING
  double t0 = GetTime();
#endif
  dcopy_(&n, dx, &incx, dy, &incy);
#ifdef FINEST_TIMING
  times[kTimeDenseFact_copy] += GetTime() - t0;
#endif
}

void callAndTime_daxpy(int n, double da, const double* dx, int incx, double* dy,
                       int incy, double* times) {
#ifdef FINEST_TIMING
  double t0 = GetTime();
#endif
  daxpy_(&n, &da, dx, &incx, dy, &incy);
#ifdef FINEST_TIMING
  times[kTimeDenseFact_axpy] += GetTime() - t0;
#endif
}

void callAndTime_dscal(int n, const double da, double* dx, int incx,
                       double* times) {
#ifdef FINEST_TIMING
  double t0 = GetTime();
#endif
  dscal_(&n, &da, dx, &incx);
#ifdef FINEST_TIMING
  times[kTimeDenseFact_scal] += GetTime() - t0;
#endif
}

int callAndTime_denseFactK(char uplo, int n, double* A, int lda,
                     const int* pivot_sign, double thresh, double* regul,
                     int* n_reg_piv, double* times) {
#ifdef FINEST_TIMING
  double t0 = GetTime();
#endif
  int info =
      denseFactK(uplo, n, A, lda, pivot_sign, thresh, regul, n_reg_piv);
#ifdef FINEST_TIMING
  times[kTimeDenseFact_fact] += GetTime() - t0;
#endif

  return info;
}