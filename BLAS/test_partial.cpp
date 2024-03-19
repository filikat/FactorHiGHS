#include <chrono>
#include <iostream>
#include <random>
#include <vector>

extern "C" void H_dpotrf(char* uplo, int* n, double* A, int* ldA, double* B,
                         int* ldB, int* info, int* k);

extern "C" int PartialFact_pos_large(int n, int k, double* A, int lda,
                                     double* B, int ldb);

extern "C" int PartialFact_ind_large(int n, int k, double* A, int lda,
                                     double* B, int ldb);

extern "C" int PartialFact_ind_small(int n, int k, double* A, int lda,
                                     double* B, int ldb);

void print(double* A, int i, int j, int m, int n, int lda) {
  for (int x = i; x < m; ++x) {
    for (int y = j; y < n; ++y) {
      printf("%10.12f ", A[x + lda * y]);
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
    A[i + i * n] += 10* n;
  }
  std::vector<double> AA(A);

  char uplo = 'l';
  int info;
  int k = n / 5;

  //print(A.data(), 0, 0, n, n, n);

  auto t0 = std::chrono::high_resolution_clock::now();
  PartialFact_ind_large(n, k, A.data(), n, &A[k + n * k], n);
  // PartialFact_pos_large(n, k, A.data(), n, &A[k + n * k], n);
  auto t1 = std::chrono::high_resolution_clock::now();
  PartialFact_ind_small(n, k, AA.data(), n, &AA[k + n * k], n);
  // H_dpotrf(&uplo, &n, AA.data(), &n, &AA[k + n * k], &n, &info, &k);
  auto t2 = std::chrono::high_resolution_clock::now();

  //print(A.data(), 0, 0, n, n, n);
  //print(AA.data(), 0, 0, n, n, n);

  double time0{};
  double time1{};
  double dur0 = std::chrono::duration<double>(t1 - t0).count();
  time0 += dur0;
  double dur1 = std::chrono::duration<double>(t2 - t1).count();
  time1 += dur1;

  double error{};
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      if (j > i) continue;
      error = std::max(error, std::fabs(A[i + j * n] - AA[i + j * n]));
    }
  }
  std::cout << "Error:" << error << '\n';
  std::cout << "Time new : " << time0 << '\n';
  std::cout << "Time old : " << time1 << '\n';

  return 0;
}