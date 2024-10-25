#include "CallAndTimeBlas.h"

#include "Blas_declaration.h"
#include "DenseFact.h"
#include "Timing.h"

void callAndTime_dsyrk(char uplo, char trans, int n, int k, double alpha,
                       const double* A, int lda, double beta, double* C,
                       int ldc, DataCollector& DC) {
#ifdef BLAS_TIMING
  double t0 = GetTime();
#endif
  dsyrk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_syrk, GetTime() - t0);
#endif
}

void callAndTime_dgemm(char transa, char transb, int m, int n, int k,
                       double alpha, const double* A, int lda, const double* B,
                       int ldb, double beta, double* C, int ldc,
                       DataCollector& DC) {
#ifdef BLAS_TIMING
  double t0 = GetTime();
#endif
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C,
         &ldc);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_gemm, GetTime() - t0);
#endif
}

void callAndTime_dtpsv(char uplo, char trans, char diag, int n,
                       const double* ap, double* x, int incx,
                       DataCollector& DC) {
#ifdef BLAS_TIMING
  double t0 = GetTime();
#endif
  dtpsv_(&uplo, &trans, &diag, &n, ap, x, &incx);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_tpsv, GetTime() - t0);
#endif
}

void callAndTime_dtrsv(char uplo, char trans, char diag, int n, const double* A,
                       int lda, double* x, int incx, DataCollector& DC) {
#ifdef BLAS_TIMING
  double t0 = GetTime();
#endif
  dtrsv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_trsv, GetTime() - t0);
#endif
}

void callAndTime_dtrsm(char side, char uplo, char trans, char diag, int m,
                       int n, double alpha, const double* A, int lda, double* B,
                       int ldb, DataCollector& DC) {
#ifdef BLAS_TIMING
  double t0 = GetTime();
#endif
  dtrsm_(&side, &uplo, &trans, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_trsm, GetTime() - t0);
#endif
}

void callAndTime_dgemv(char trans, int m, int n, double alpha, const double* A,
                       int lda, const double* x, int incx, double beta,
                       double* y, int incy, DataCollector& DC) {
#ifdef BLAS_TIMING
  double t0 = GetTime();
#endif
  dgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_gemv, GetTime() - t0);
#endif
}

void callAndTime_dcopy(int n, const double* dx, int incx, double* dy, int incy,
                       DataCollector& DC) {
#ifdef BLAS_TIMING
  double t0 = GetTime();
#endif
  dcopy_(&n, dx, &incx, dy, &incy);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_copy, GetTime() - t0);
#endif
}

void callAndTime_daxpy(int n, double da, const double* dx, int incx, double* dy,
                       int incy, DataCollector& DC) {
#ifdef BLAS_TIMING
  double t0 = GetTime();
#endif
  daxpy_(&n, &da, dx, &incx, dy, &incy);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_axpy, GetTime() - t0);
#endif
}

void callAndTime_dscal(int n, const double da, double* dx, int incx,
                       DataCollector& DC) {
#ifdef BLAS_TIMING
  double t0 = GetTime();
#endif
  dscal_(&n, &da, dx, &incx);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_scal, GetTime() - t0);
#endif
}

int callAndTime_denseFactK(char uplo, int n, double* A, int lda,
                           const int* pivot_sign, double thresh, double* regul,
                           DataCollector& DC) {
#ifdef FINE_TIMING
  double t0 = GetTime();
#endif
  int info = denseFactK(uplo, n, A, lda, pivot_sign, thresh, regul, DC);
#ifdef FINE_TIMING
  DC.sumTime(kTimeDenseFact_fact, GetTime() - t0);
#endif

  return info;
}