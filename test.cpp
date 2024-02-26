#include <fstream>
#include <iostream>

#include "Analyze.h"

int main() {
  std::ifstream file("matlab/data.txt");

  int n{};
  file >> n;

  int nz{};
  file >> nz;

  std::vector<int> ptr(n + 1);
  for (int i = 0; i < n + 1; ++i) {
    file >> ptr[i];
  }

  std::vector<int> rows(nz);
  for (int i = 0; i < nz; ++i) {
    file >> rows[i];
  }

  file.close();

  Analyze An(rows.data(), ptr.data(), n, nz, false);
  An.GetPermutation();

  for (int i = 0; i < n; ++i) {
    std::cout << An.metis_perm[i] << ' ';
  }
  std::cout << "\n\n";

  return 0;
}
