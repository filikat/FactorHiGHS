#include <chrono>
#include <iostream>
#include <random>
#include <vector>

extern "C" void dgemm(char* transa, char* transb, int* m, int* n, int* k,
                      double* alpha, double* A, int* lda, double* B, int* ldb,
                      double* beta, double* C, int* ldc);

extern "C" void H_dgemm_block(char* transa, char* transb, int* m, int* n,
                              int* k, double* alpha, double* A, int* lda,
                              double* B, int* ldb, double* beta, double* C,
                              int* ldc);

void test(char transa, char transb, int m, int n, int k, double alpha,
          double* A, int lda, double* B, int ldb, double beta, double* C,
          double* HC, int ldc, double& time0, double& time1) {
  auto t0 = std::chrono::high_resolution_clock::now();
  dgemm(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
  auto t1 = std::chrono::high_resolution_clock::now();
  H_dgemm_block(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
                HC, &ldc);
  auto t2 = std::chrono::high_resolution_clock::now();

  double dur0 = std::chrono::duration<double>(t1 - t0).count();
  time0 += dur0;
  double dur1 = std::chrono::duration<double>(t2 - t1).count();
  time1 += dur1;

  double error = 0.0;
  double norm = 0.0;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      error = std::max(error, std::fabs(HC[j * m + i] - C[j * m + i]));
      norm += std::fabs(C[j * m + i]);
    }
  }
  error /= norm;

  printf("Test %c%c: error %8.1e, B %5.3f, H %5.3f, ratio %5.1f", transa,
         transb, error, dur0, dur1, dur1 / dur0);
  if (error > 1e-6) {
    printf("FAILED <===\n");
  } else {
    printf("\n");
  }
}

int main() {
  int m = 2000;
  int n = 2000;
  int k = 1000;

  double time0{};
  double time1{};

  std::random_device rand_dev;
  std::mt19937 rng(rand_dev());
  std::uniform_real_distribution<double> distr(-50.0, 50.0);

  std::vector<double> A(m * k);
  for (int i = 0; i < A.size(); ++i) {
    A[i] = distr(rng);
  }

  std::vector<double> B(m * k);
  for (int i = 0; i < B.size(); ++i) {
    B[i] = distr(rng);
  }

  std::vector<double> C(m * n);
  for (int i = 0; i < C.size(); ++i) {
    C[i] = distr(rng);
  }
  std::vector<double> HC(C);

  double alpha = distr(rng);
  double beta = distr(rng);

  std::vector<double> CC;
  std::vector<double> HCC;

  CC = C;
  HCC = HC;
  test('N', 'N', m, n, k, alpha, A.data(), m, B.data(), k, beta, CC.data(),
       HCC.data(), m, time0, time1);

  CC = C;
  HCC = HC;
  test('T', 'N', m, n, k, alpha, A.data(), k, B.data(), k, beta, CC.data(),
       HCC.data(), m, time0, time1);

  CC = C;
  HCC = HC;
  test('N', 'T', m, n, k, alpha, A.data(), m, B.data(), n, beta, CC.data(),
       HCC.data(), m, time0, time1);

  CC = C;
  HCC = HC;
  test('T', 'T', m, n, k, alpha, A.data(), k, B.data(), n, beta, CC.data(),
       HCC.data(), m, time0, time1);

  printf("Time BLAS %5.2f\n", time0);
  printf("Time H %5.2f\n", time1);
  printf("Ratio %5.2f\n", time1 / time0);

  return 0;
}