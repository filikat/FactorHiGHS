#include <chrono>
#include <iostream>
#include <random>
#include <vector>

extern "C" void dpotrf(char* uplo, int* n, double* A, int* ldA, int* info);
extern "C" void H_dpotrf(char* uplo, int* n, double* A, int* ldA, int* info);
extern "C" void H_dpotrf_packed(int* N, double* A, int* NB);

void pack(const std::vector<double>& A, std::vector<double>& P, int n, int nb);
void unpack(const std::vector<double>& P, std::vector<double>& A, int n,
            int nb);

void test(char uplo, int n, const std::vector<double>& A, int lda) {
  std::vector<double> AA(A);
  std::vector<double> HA(A);
  std::vector<double> PA(A);

  double time_lapack{};
  double time_packing{};
  double time_H{};
  double time_H_packed{};

  auto t0 = std::chrono::high_resolution_clock::now();
  int info{};
  dpotrf(&uplo, &n, AA.data(), &lda, &info);
  auto t1 = std::chrono::high_resolution_clock::now();
  if (info != 0) {
    std::cout << "Wrong info\n";
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  H_dpotrf(&uplo, &n, HA.data(), &lda, &info);
  auto t3 = std::chrono::high_resolution_clock::now();
  if (info != 0) {
    std::cout << "Wrong info\n";
  }

  int nb = 128;
  std::vector<double> P(n * (n + 1) / 2);
  auto t4 = std::chrono::high_resolution_clock::now();
  pack(PA, P, n, nb);
  auto t5 = std::chrono::high_resolution_clock::now();
  H_dpotrf_packed(&n, P.data(), &nb);
  auto t6 = std::chrono::high_resolution_clock::now();
  unpack(P, PA, n, nb);
  auto t7 = std::chrono::high_resolution_clock::now();

  double dur0 = std::chrono::duration<double>(t1 - t0).count();
  time_lapack += dur0;
  double dur1 = std::chrono::duration<double>(t3 - t2).count();
  time_H += dur1;
  double dur2 = std::chrono::duration<double>(t5 - t4).count();
  time_packing += dur2;
  double dur3 = std::chrono::duration<double>(t6 - t5).count();
  time_H_packed += dur3;
  double dur4 = std::chrono::duration<double>(t7 - t6).count();
  time_packing += dur4;

  double error = 0.0;
  double norm = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      error = std::max(error, std::fabs(HA[j * n + i] - AA[j * n + i]));
      error = std::max(error, std::fabs(PA[j * n + i] - AA[j * n + i]));
      norm += std::fabs(AA[j * n + i]);
    }
  }
  error /= norm;

  printf("Error %8.1e\n", error);
  printf("Times:\n");
  printf("%10s %5.3f\n", "Lapack", time_lapack);
  printf("%10s %5.3f\n", "H", time_H);
  printf("%10s %5.3f\n", "H packed", time_H_packed);
  printf("%10s %5.3f\n", "packing", time_packing);
  printf("%10s %5.3f\n", "packed tot", time_H_packed + time_packing);
  if (error > 1e-6) {
    printf("FAILED <===\n");
  }
}

int main() {
  int n = 30000;

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

  test('L', n, A, n);

  return 0;
}