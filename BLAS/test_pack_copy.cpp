#include <stdlib.h>

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "../DenseFact_declaration.h"

std::random_device rand_dev;
std::mt19937 generator(rand_dev());
std::uniform_int_distribution<int> distr(1, 20);

int main() {
  // ===========================================================================
  // Set up
  // ===========================================================================
  int nb = 5;
  int nrow = 18;
  int ncol = 7;

  int l = nrow * ncol - ncol * (ncol - 1) / 2;

  std::vector<double> A(l);

  int next{};
  for (int j = 0; j < ncol; ++j) {
    for (int i = j; i < nrow; ++i) {
      A[next] = distr(generator);
      if (i == j) A[next] += nrow * 20;
      ++next;
    }
  }

  // ===========================================================================
  // Create and print full A
  // ===========================================================================

  if (nrow < 20) {
    std::vector<double> fullA(nrow * ncol);
    int startA{};
    for (int j = 0; j < ncol; ++j) {
      for (int i = j; i < nrow; ++i) {
        fullA[i + nrow * j] = A[startA++];
      }
    }

    for (int i = 0; i < nrow; ++i) {
      for (int j = 0; j < ncol; ++j) {
        printf("%5.0f ", fullA[i + nrow * j]);
      }
      printf("\n");
    }
  }

  auto t0 = std::chrono::high_resolution_clock::now();
  PackedToHybrid(A.data(), nrow, ncol, nb);
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> d01 = t1 - t0;
  printf("%f\n\n", d01.count());

  // for (double d : A) printf("%.0f ", d);
  // printf("\n");

  int lb = nrow - ncol;
  lb = lb * (lb + 1) / 2;
  std::vector<double> B(lb);

  std::vector<double> times(TimesInd::t_size);
  PartialFactPosPacked(nrow, ncol, A.data(), nb, B.data(), times.data());


  for (double d : B) {
    printf("%f ", d);
  }
  printf("\n");

  /*
    //
    ===========================================================================
    // Buffer
    //
    ===========================================================================
    std::vector<double> buf(n * nb);

    // copy A into buf, with spaces
    int start{};
    int pos{};
    for (int j = 0; j < nb; ++j) {
      int N = n - j;
      int i_one = 1;
      dcopy(&N, &A[start], &i_one, &buf[pos], &i_one);
      start += N;
      pos += n + 1;
    }

    //
    ===========================================================================
    // Method 1
    //
    ===========================================================================
    std::vector<double> P1(l);

    auto t0 = std::chrono::high_resolution_clock::now();
    // copy buf into P1 by rows
    start = 0;
    for (int i = 0; i < n; ++i) {
      int N = std::min(i + 1, nb);
      int i_one = 1;
      dcopy(&N, &buf[i], &n, &P1[start], &i_one);
      start += N;
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    //
    ===========================================================================
    // Method 2
    //
    ===========================================================================
    std::vector<double> P2(l);

    auto t2 = std::chrono::high_resolution_clock::now();
    int start_buf{};
    int start_P_orig{};
    for (int j = 0; j < nb; ++j) {
      int jump = j + 1;
      int start_P = start_P_orig;
      start_P_orig += j + 2;

      while (jump < nb) {
        P2[start_P] = buf[start_buf];
        ++start_buf;
        start_P += jump;
        ++jump;
      }

      int N = n - nb + 1;
      int i_one = 1;
      dcopy(&N, &buf[start_buf], &i_one, &P2[start_P], &nb);
      start_buf += N + j + 1;
    }
    auto t3 = std::chrono::high_resolution_clock::now();

    //
    ===========================================================================
    // Method 3
    //
    ===========================================================================
    auto t4 = std::chrono::high_resolution_clock::now();
    std::vector<double> P3(l);
    pack(buf, P3, n, nb, nb);
    auto t5 = std::chrono::high_resolution_clock::now();

    //
    ===========================================================================
    // Results
    //
    ===========================================================================
    for (int i = 0; i < l; ++i) {
      if (P1[i] != P2[i] || P1[i] != P3[i]) {
        printf("Error\n");
      }
    }

    std::chrono::duration<double> d01 = t1 - t0;
    std::chrono::duration<double> d23 = t3 - t2;
    std::chrono::duration<double> d45 = t5 - t4;
    printf("Method 1: %f\n", d01.count());
    printf("Method 2: %f\n", d23.count());
    printf("Method 3: %f\n", d45.count());
  */
  return 0;
}
