#include "Symbolic.h"

#include <iostream>

void Symbolic::Print() const {
  printf("Symbolic factorisation:\n");
  printf(" - type                 %s\n",
         type == FactType::NormEq ? "Normal equations" : "Augmented system");
  printf(" - size                 %d\n", n);
  printf(" - nonzero entries      %.2e\n", (double)nz);
  printf(" - density              %.2f\n", ((double)nz / n) / n);
  printf(" - fill in              %.2f\n", fillin);
  printf(" - supernodes           %d\n", sn);
  printf(" - largest supernode    %d\n", largestSn);
  printf(" - largest front        %d\n", largestFront);
  printf(" - dense operations     %.2e\n", operations);
  printf(" - assembly operations  %.2e\n", assemblyOp);
  printf(" - artificial nonzeros  %.2e (%4.1f%%)\n", (double)artificialNz,
         (double)artificialNz / nz * 100);
  printf(" - artificial ops       %.2e (%4.1f%%)\n", artificialOp,
         artificialOp / operations * 100);
}

FactType Symbolic::Type() const { return type; }
int Symbolic::Size() const { return n; }
int Symbolic::Nz() const { return nz; }
double Symbolic::Ops() const { return operations; }
double Symbolic::AssemblyOps() const { return assemblyOp; }
int Symbolic::Sn() const { return sn; }
int Symbolic::Rows(int i) const { return rows[i]; }
int Symbolic::Ptr(int i) const { return ptr[i]; }
int Symbolic::SnStart(int i) const { return snStart[i]; }
int Symbolic::RelindCols(int i) const { return relindCols[i]; }
int Symbolic::RelindClique(int i, int j) const { return relindClique[i][j]; }
int Symbolic::ConsecutiveSums(int i, int j) const {
  return consecutiveSums[i][j];
}

const std::vector<int>& Symbolic::Ptr() const { return ptr; }
const std::vector<int>& Symbolic::Perm() const { return perm; }
const std::vector<int>& Symbolic::Iperm() const { return iperm; }
const std::vector<int>& Symbolic::SnParent() const { return snParent; }
const std::vector<int>& Symbolic::SnStart() const { return snStart; }
