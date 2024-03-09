#ifndef HBLAS_H
#define HBLAS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define max(i, j) ((i) >= (j) ? (i) : (j))

void H_dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha, double* a,
             int* lda, double* beta, double* c, int* ldc);

void H_dtrsm(char* side, char* uplo, char* trans, char* diag, int* m, int* n,
             double* alpha, double* a, int* lda, double* b, int* ldb);

void H_dgemm(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
             double* A, int* lda, double* B, int* ldb, double* beta, double* C,
             int* ldc);

void H_dgemm_block(char* transa, char* transb, int* m, int* n, int* k,
                   double* alpha, double* A, int* lda, double* B, int* ldb,
                   double* beta, double* C, int* ldc);

void H_dgepp(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
             double* A, int* ldA, double* B, int* ldB, double* C, int* ldC);

void H_dgebp(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
             double* A, int* ldA, double* B, double* C, int* ldC);

void H_kernel(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
              double* A, double* B, double* C, int* ldC);

// pagesize (16MB) / 16
const int kc = 512;

// min(L2 size, TLB size) / (16 * kc) = min(4MB, 32MB) / (16 * 1024) ??
// Unsure about L2 and TLB
const int mc = 256;

// small, depending on flops/cycle and bandwidth of cache
const int mr = 16;
const int nr = 16;

double time_copy_B = 0.0;
double time_copy_A = 0.0;
double time_kernel = 0.0;

double get_time() {
  struct timespec now;
  clock_gettime(CLOCK_REALTIME, &now);
  return now.tv_sec + now.tv_nsec * 1e-9;
}

#endif