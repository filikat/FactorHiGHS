#include <chrono>
#include <iostream>
#include <random>
#include <vector>

extern "C" void dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha,
                      double* a, int* lda, double* beta, double* c, int* ldc);

extern "C" void H_dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha,
                        double* a, int* lda, double* beta, double* c, int* ldc);

void test(char uplo, char trans, int n, int k, double alpha, double* A, int lda,
          double beta, double* C, double* HC, int ldc, double& time0,
          double& time1) {
  auto t0 = std::chrono::high_resolution_clock::now();
  dsyrk(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc);
  auto t1 = std::chrono::high_resolution_clock::now();
  H_dsyrk(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, HC, &ldc);
  auto t2 = std::chrono::high_resolution_clock::now();

  double dur0 = std::chrono::duration<double>(t1 - t0).count();
  time0 += dur0;
  double dur1 = std::chrono::duration<double>(t2 - t1).count();
  time1 += dur1;

  double error = 0.0;
  double norm = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      error = std::max(error, std::fabs(HC[j * n + i] - C[j * n + i]));
      norm += std::fabs(C[j * n + i]);
    }
  }
  error /= norm;

  printf("Test %c%c: error %8.1e, B %5.3f, H %5.3f, ratio %5.1f", uplo, trans,
         error, dur0, dur1, dur1 / dur0);
  if (error > 1e-6) {
    printf("FAILED <===\n");
  } else {
    printf("\n");
  }
}

int main() {
  int n = 10000;
  int k = 1000;

  double time0{};
  double time1{};

  std::random_device rand_dev;
  std::mt19937 rng(rand_dev());
  std::uniform_real_distribution<double> distr(-50.0, 50.0);

  std::vector<double> A(n * k);
  for (int i = 0; i < k * n; ++i) {
    A[i] = distr(rng);
  }

  double alpha = distr(rng);
  double beta = distr(rng);

  std::vector<double> C(n * n);
  std::vector<double> HC(n * n);

  test('U', 'N', n, k, alpha, A.data(), n, beta, C.data(), HC.data(), n, time0,
       time1);
  test('L', 'N', n, k, alpha, A.data(), n, beta, C.data(), HC.data(), n, time0,
       time1);

  C.resize(k * k);
  HC.resize(k * k);
  test('U', 'T', k, n, alpha, A.data(), n, beta, C.data(), HC.data(), k, time0,
       time1);
  test('L', 'T', k, n, alpha, A.data(), n, beta, C.data(), HC.data(), k, time0,
       time1);

  printf("Time BLAS %5.2f\n", time0);
  printf("Time H %5.2f\n", time1);
  printf("Ratio %5.2f\n", time1 / time0);

  return 0;
}