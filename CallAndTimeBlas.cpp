#include "CallAndTimeBlas.h"

#include "Auxiliary.h"
#include "Blas_declaration.h"
#include "DenseFact.h"
#include "Timing.h"

void callAndTime_dsyrk(char uplo, char trans, int n, int k, double alpha,
                       const double* A, int lda, double beta, double* C,
                       int ldc, DataCollector& DC) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  dsyrk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_syrk, clock.stop());
#endif
}

void callAndTime_dgemm(char transa, char transb, int m, int n, int k,
                       double alpha, const double* A, int lda, const double* B,
                       int ldb, double beta, double* C, int ldc,
                       DataCollector& DC) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C,
         &ldc);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_gemm, clock.stop());
#endif
}

void callAndTime_dtpsv(char uplo, char trans, char diag, int n,
                       const double* ap, double* x, int incx,
                       DataCollector& DC) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  dtpsv_(&uplo, &trans, &diag, &n, ap, x, &incx);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_tpsv, clock.stop());
#endif
}

void callAndTime_dtrsv(char uplo, char trans, char diag, int n, const double* A,
                       int lda, double* x, int incx, DataCollector& DC) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  dtrsv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_trsv, clock.stop());
#endif
}

void callAndTime_dtrsm(char side, char uplo, char trans, char diag, int m,
                       int n, double alpha, const double* A, int lda, double* B,
                       int ldb, DataCollector& DC) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  dtrsm_(&side, &uplo, &trans, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_trsm, clock.stop());
#endif
}

void callAndTime_dgemv(char trans, int m, int n, double alpha, const double* A,
                       int lda, const double* x, int incx, double beta,
                       double* y, int incy, DataCollector& DC) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  dgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_gemv, clock.stop());
#endif
}

void callAndTime_dger(int m, int n, double alpha, const double* x, int incx,
                      const double* y, int incy, double* A, int lda,
                      DataCollector& DC) {
  dger_(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
}

void callAndTime_dcopy(int n, const double* dx, int incx, double* dy, int incy,
                       DataCollector& DC) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  dcopy_(&n, dx, &incx, dy, &incy);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_copy, clock.stop());
#endif
}

void callAndTime_daxpy(int n, double da, const double* dx, int incx, double* dy,
                       int incy, DataCollector& DC) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  daxpy_(&n, &da, dx, &incx, dy, &incy);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_axpy, clock.stop());
#endif
}

void callAndTime_dscal(int n, const double da, double* dx, int incx,
                       DataCollector& DC) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  dscal_(&n, &da, dx, &incx);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_scal, clock.stop());
#endif
}

void callAndTime_dswap(int n, double* dx, int incx, double* dy, int incy,
                       DataCollector& DC) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  dswap_(&n, dx, &incx, dy, &incy);
#ifdef BLAS_TIMING
  DC.sumTime(kTimeBlas_scal, clock.stop());
#endif
}

int callAndTime_denseFactK(char uplo, int n, double* A, int lda,
                           const int* pivot_sign, double thresh, double* regul,
                           DataCollector& DC, int sn, int bl) {
#ifdef FINE_TIMING
  Clock clock;
#endif
  int info = denseFactK(uplo, n, A, lda, pivot_sign, thresh, regul, DC, sn, bl);
#ifdef FINE_TIMING
  DC.sumTime(kTimeDenseFact_fact, clock.stop());
#endif

  return info;
}