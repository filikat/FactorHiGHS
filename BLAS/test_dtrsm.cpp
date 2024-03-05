#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "Hblas.h"

extern "C" void dtrsm(char* side, char* uplo, char* trans, char* diag, int* m,
                      int* n, double* alpha, double* A, int* lda, double* B,
                      int* ldb);

void test(char side, char uplo, char trans, char diag, int m, int n,
          double alpha, double* A, int lda, double* B, double* HB, int ldb,
          double& time0, double& time1) {
  auto t0 = std::chrono::high_resolution_clock::now();
  dtrsm(&side, &uplo, &trans, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
  auto t1 = std::chrono::high_resolution_clock::now();
  H_dtrsm(&side, &uplo, &trans, &diag, &m, &n, &alpha, A, &lda, HB, &ldb);
  auto t2 = std::chrono::high_resolution_clock::now();

  double dur0 = std::chrono::duration<double>(t1 - t0).count();
  time0 += dur0;
  double dur1 = std::chrono::duration<double>(t2 - t1).count();
  time1 += dur1;

  double error = 0.0;
  double norm = 0.0;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      error = std::max(error, std::fabs(HB[j * m + i] - B[j * m + i]));
      norm += std::fabs(B[j * m + i]);
    }
  }
  error /= norm;

  printf("Test %c%c%c%c: error %8.1e, B %5.3f, H %5.3f, ratio %5.1f", side,
         uplo, trans, diag, error, dur0, dur1, dur1 / dur0);
  if (error > 1e-6) {
    printf("FAILED <===\n");
  } else {
    printf("\n");
  }
}

int main() {
  int n = 10000;
  int m = 100;

  double time0{};
  double time1{};

  std::random_device rand_dev;
  std::mt19937 rng(rand_dev());
  std::uniform_real_distribution<double> distr(-50.0, 50.0);

  std::vector<double> B(m * n);
  std::vector<double> HB(m * n);
  for (int i = 0; i < B.size(); ++i) {
    B[i] = distr(rng);
    HB[i] = B[i];
  }

  std::vector<double> BB;
  std::vector<double> HBB;

  std::vector<double> A(m * m);
  for (int i = 0; i < A.size(); ++i) {
    A[i] = distr(rng);
  }

  double alpha = distr(rng);

  char side = 'L';
  char uplo[2] = {'U', 'L'};
  char trans[2] = {'N', 'T'};
  char diag[2] = {'N', 'U'};

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        BB = B;
        HBB = HB;
        test(side, uplo[i], trans[j], diag[k], m, n, alpha, A.data(), m,
             BB.data(), HBB.data(), m, time0, time1);
      }
    }
  }

  side = 'R';
  A.resize(n * n);
  for (int i = 0; i < A.size(); ++i) {
    A[i] = distr(rng);
  }

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        BB = B;
        HBB = HB;
        test(side, uplo[i], trans[j], diag[k], m, n, alpha, A.data(), n,
             BB.data(), HBB.data(), m, time0, time1);
      }
    }
  }

  printf("Time BLAS %5.2f\n", time0);
  printf("Time H %5.2f\n", time1);
  printf("Ratio %5.2f\n", time1 / time0);

  return 0;
}