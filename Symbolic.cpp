#include "Symbolic.h"

#include <iostream>

void Symbolic::Print() const {
  printf("Symbolic factorisation:\n");
  printf(" - size                 %d\n", n);
  printf(" - nonzero entries      %.2e\n", (double)nz);
  printf(" - fill in              %.2f\n", fillin);
  printf(" - supernodes           %d\n", sn);
  printf(" - largest supernode    %d\n", largestSn);
  printf(" - largest front        %d\n", largestFront);
  printf(" - dense operations     %.2e\n", operations);
  printf(" - assembly operations  %.2e\n", assemblyOp);
  printf(" - artificial nonzeros  %.2e (%2.1f%%)\n", (double)artificialNz,
         (double)artificialNz / nz * 100);
  printf(" - artificial ops       %.2e (%2.1f%%)\n", artificialOp,
         artificialOp / operations * 100);
}

int Symbolic::Size() const { return n; }
int Symbolic::Nz() const { return nz; }
double Symbolic::Ops() const { return operations; }
double Symbolic::AssemblyOps() const { return assemblyOp; }
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
