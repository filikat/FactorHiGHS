#include <chrono>
#include <iostream>
#include <random>
#include <vector>

extern "C" void dpotrf(char* uplo, int* n, double* A, int* ldA, int* info);

extern "C" void H_dpotrf(char* uplo, int* n, double* A, int* ldA, int* info);

void test(char uplo, int n, const std::vector<double>& A, int lda,
          double& time0, double& time1) {
  std::vector<double> AA(A);
  std::vector<double> HA(A);

  auto t0 = std::chrono::high_resolution_clock::now();
  int info{};
  dpotrf(&uplo, &n, AA.data(), &lda, &info);
  if (info != 0) {
    std::cout << "Wrong info\n";
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  H_dpotrf(&uplo, &n, HA.data(), &lda, &info);
  auto t2 = std::chrono::high_resolution_clock::now();

  double dur0 = std::chrono::duration<double>(t1 - t0).count();
  time0 += dur0;
  double dur1 = std::chrono::duration<double>(t2 - t1).count();
  time1 += dur1;

  double error = 0.0;
  double norm = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      error = std::max(error, std::fabs(HA[j * n + i] - AA[j * n + i]));
      norm += std::fabs(AA[j * n + i]);
    }
  }
  error /= norm;

  printf("Test %c: error %8.1e, B %5.3f, H %5.3f, ratio %5.1f", uplo, error,
         dur0, dur1, dur1 / dur0);
  if (error > 1e-6) {
    printf("FAILED <===\n");
  } else {
    printf("\n");
  }
}

int main() {
  int n = 10000;

  double time0{};
  double time1{};

  std::random_device rand_dev;
  std::mt19937 rng(rand_dev());
  std::uniform_real_distribution<double> distr(-50.0, 50.0);

  std::vector<double> A(n * n);
  for (int i = 0; i < A.size(); ++i) {
    A[i] = distr(rng);
  }

  for (int i = 0; i < n; ++i) {
    A[i + i * n] += 50 * n;
  }

  test('L', n, A, n, time0, time1);
  test('U', n, A, n, time0, time1);

  printf("Time BLAS %5.2f\n", time0);
  printf("Time H %5.2f\n", time1);
  printf("Ratio %5.2f\n", time1 / time0);

  return 0;
}