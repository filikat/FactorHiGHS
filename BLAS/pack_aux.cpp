#include <iostream>
#include <vector>

void pack(const std::vector<double>& A, std::vector<double>& P, int n, int nb) {
  int start0{};
  int start{};

  for (int k = 0; k <= (n - 1) / nb; ++k) {
    int block_size = std::min(nb, n - k * nb);

    for (int j = 0; j < block_size; ++j) {
      start = start0 + (j + 1) * (j + 2) / 2 - 1;
      int jump = j + 1;
      int coloffset = (j + k * nb) * n;

      for (int i = k * nb + j; i < n; ++i) {
        P[start] = A[i + coloffset];
        start += jump;
        jump = jump < nb ? jump + 1 : nb;
      }
    }

    start0 = start - nb + 1;
  }
}

void unpack(const std::vector<double>& P, std::vector<double>& A, int n,
            int nb) {
  int start0{};
  int start{};

  for (int k = 0; k <= (n - 1) / nb; ++k) {
    int block_size = std::min(nb, n - k * nb);

    for (int j = 0; j < block_size; ++j) {
      start = start0 + (j + 1) * (j + 2) / 2 - 1;
      int jump = j + 1;

      for (int i = k * nb + j; i < n; ++i) {
        A[i + (j + k * nb) * n] = P[start];
        start += jump;
        jump = jump < nb ? jump + 1 : nb;
      }
    }
    start0 = start - nb + 1;
  }
}