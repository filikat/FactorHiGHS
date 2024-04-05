#include "Symbolic.h"

#include <iostream>

void Symbolic::Print() const {
  printf("Symbolic factorization:\n");
  printf(" - size                %d\n", n);
  printf(" - nonzero entries     %.2e\n", (double)nz);
  printf(" - supernodes found    %d\n", sn);
  printf(" - artificial nonzeros %d (%2.1f%%)\n", artificialNz,
         (double)artificialNz / nz * 100);
  printf(" - operations count    %.2e\n", operations);
  printf(" - largest front       %d\n", largestFront);
  printf(" - largest supernode   %d\n", largestSn);
}

int Symbolic::Size() const { return n; }
int Symbolic::Nz() const { return nz; }
int Symbolic::Ops() const { return operations; }
int Symbolic::Sn() const { return sn; }
int Symbolic::Rows(int i) const { return rows[i]; }
int Symbolic::Ptr(int i) const { return ptr[i]; }
int Symbolic::Sn_start(int i) const { return sn_start[i]; }
int Symbolic::Relind_cols(int i) const { return relind_cols[i]; }
int Symbolic::Relind_clique(int i, int j) const { return relind_clique[i][j]; }
int Symbolic::ConsecutiveSums(int i, int j) const {
  return consecutiveSums[i][j];
}

const std::vector<int>& Symbolic::Ptr() const { return ptr; }
const std::vector<int>& Symbolic::Perm() const { return perm; }
const std::vector<int>& Symbolic::Iperm() const { return iperm; }
const std::vector<int>& Symbolic::Sn_parent() const { return sn_parent; }
const std::vector<int>& Symbolic::Sn_start() const { return sn_start; }
