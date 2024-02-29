#include <fstream>
#include <iostream>

#include "Analyze.h"

int main() {
  // Read the problem
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

  // Symbolic factorization
  Symbolic S;

  Analyze An(rows.data(), ptr.data(), n, nz, true);
  An.Run(S, 1, 1);

  // Write to file

  std::ofstream out_file;

  out_file.open("matlab/perm.txt");
  for (int i : S.perm) {
    out_file << i << '\n';
  }
  out_file.close();

  out_file.open("matlab/parent.txt");
  for (int i : S.parent) {
    out_file << i << '\n';
  }
  out_file.close();

  out_file.open("matlab/rowcount.txt");
  for (int i : S.rowcount) {
    out_file << i << '\n';
  }
  out_file.close();

  out_file.open("matlab/colcount.txt");
  for (int i : S.colcount) {
    out_file << i << '\n';
  }
  out_file.close();

  out_file.open("matlab/ptrL.txt");
  for (int i : S.ptr) {
    out_file << i << '\n';
  }
  out_file.close();

  out_file.open("matlab/rowsL.txt");
  for (int i : S.rows) {
    out_file << i << '\n';
  }
  out_file.close();

  out_file.open("matlab/fsn_ptr.txt");
  for (int i : S.fsn_ptr) {
    out_file << i << '\n';
  }
  out_file.close();

  return 0;
}
