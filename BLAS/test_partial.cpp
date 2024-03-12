#include <chrono>
#include <iostream>
#include <random>
#include <vector>

extern "C" void H_dpotrf(char* uplo, int* n, double* A, int* ldA, int* info,
                         int* k);
extern "C" void H_dpotrf_right(char* uplo, int* n, double* A, int* ldA,
                               int* info, int* k);

void print(double* A, int i, int j, int m, int n, int lda) {
  for (int x = i; x < m; ++x) {
    for (int y = j; y < n; ++y) {
      printf("%9.3f ", A[x + lda * y]);
    }
    printf("\n");
  }
  printf("\n\n");
}

int main() {
  int n = 10000;

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
  std::vector<double> AA(A);

  char uplo = 'U';
  int info;
  int k = n;

  auto t0 = std::chrono::high_resolution_clock::now();
  H_dpotrf(&uplo, &n, A.data(), &n, &info, &k);
  auto t1 = std::chrono::high_resolution_clock::now();
  H_dpotrf_right(&uplo, &n, AA.data(), &n, &info, &k);
  auto t2 = std::chrono::high_resolution_clock::now();

  double time0{};
  double time1{};
  double dur0 = std::chrono::duration<double>(t1 - t0).count();
  time0 += dur0;
  double dur1 = std::chrono::duration<double>(t2 - t1).count();
  time1 += dur1;

  double error{};
  for (int i = 0; i < A.size(); ++i) {
    error = std::max(error, std::fabs(A[i] - AA[i]));
  }
  std::cout << "Error:" << error << '\n';
  std::cout << "Time left : " << time0 << '\n';
  std::cout << "Time right: " << time1 << '\n';

  return 0;
}