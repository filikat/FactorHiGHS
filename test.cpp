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

  Analyze An(rows.data(), ptr.data(), n, nz, true);
  An.GetPermutation();
  An.Permute(An.iperm);
  An.ETree();
  An.Postorder();
  An.RowColCount();

  std::vector<int> rowsL(An.nzL);
  std::vector<int> ptrL(n + 1);
  An.ColPattern(rowsL, ptrL);

  std::ofstream out_file;

  out_file.open("matlab/perm.txt");
  for (int i = 0; i < n; ++i) {
    out_file << An.perm[i] << '\n';
  }
  out_file.close();

  out_file.open("matlab/parent.txt");
  for (int i = 0; i < n; ++i) {
    out_file << An.parent[i] << '\n';
  }
  out_file.close();

  out_file.open("matlab/rowcount.txt");
  for (int i = 0; i < n; ++i) {
    out_file << An.rowcount[i] << '\n';
  }
  out_file.close();

  out_file.open("matlab/colcount.txt");
  for (int i = 0; i < n; ++i) {
    out_file << An.colcount[i] << '\n';
  }
  out_file.close();

  out_file.open("matlab/ptrL.txt");
  for (int i = 0; i < n + 1; ++i) {
    out_file << ptrL[i] << '\n';
  }
  out_file.close();

  out_file.open("matlab/rowsL.txt");
  for (int i = 0; i < An.nzL; ++i) {
    out_file << rowsL[i] << '\n';
  }
  out_file.close();

  return 0;
}
