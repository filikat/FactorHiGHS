#include "Symbolic.h"

#include <iostream>

Symbolic::Symbolic(FormatType format_type) : format_type_{format_type} {}

FormatType Symbolic::formatType() const { return format_type_; }
int Symbolic::blockSize() const { return block_size_; }
int Symbolic::size() const { return n_; }
int Symbolic::nz() const { return nz_; }
double Symbolic::ops() const { return dense_ops_; }
double Symbolic::assemblyOps() const { return assembly_ops_; }
int Symbolic::sn() const { return sn_; }
int Symbolic::rows(int i) const { return rows_[i]; }
int Symbolic::ptr(int i) const { return ptr_[i]; }
int Symbolic::snStart(int i) const { return sn_start_[i]; }
int Symbolic::relindCols(int i) const { return relind_cols_[i]; }
int Symbolic::relindClique(int i, int j) const { return relind_clique_[i][j]; }
int Symbolic::consecutiveSums(int i, int j) const {
  return consecutive_sums_[i][j];
}
int Symbolic::stackSize() const { return max_stack_entries_; }
int Symbolic::maxCliqueSize() const { return max_clique_entries_; }

const std::vector<int>& Symbolic::ptr() const { return ptr_; }
const std::vector<int>& Symbolic::iperm() const { return iperm_; }
const std::vector<int>& Symbolic::snParent() const { return sn_parent_; }
const std::vector<int>& Symbolic::snStart() const { return sn_start_; }
const std::vector<int>& Symbolic::pivotSign() const { return pivot_sign_; }

void printMemory(double mem) {
  if (mem < 1024)
    printf("%.1f B\n", mem);
  else if (mem < 1024 * 1024)
    printf("%.1f KB\n", mem / 1024);
  else if (mem < 1024 * 1024 * 1024)
    printf("%.1f MB\n", mem / 1024 / 1024);
  else
    printf("%.1f GB\n", mem / 1024 / 1024 / 1024);
}

void Symbolic::print() const {
  printf("Symbolic factorisation:\n");
  printf(" - size                 %d\n", n_);
  printf(" - nonzero entries      %.2e\n", nz_);
  printf(" - density              %.2f\n", (nz_ / n_) / n_);
  printf(" - fill in              %.2f\n", fillin_);
  printf(" - supernodes           %d\n", sn_);
  printf(" - largest supernode    %d\n", largest_sn_);
  printf(" - largest front        %d\n", largest_front_);
  printf(" - dense operations     %.2e\n", dense_ops_);
  printf(" - assembly operations  %.2e\n", assembly_ops_);
  printf(" - artificial nonzeros  %.2e (%4.1f%%)\n", (double)artificial_nz_,
         (double)artificial_nz_ / nz_ * 100);
  printf(" - artificial ops       %.2e (%4.1f%%)\n", artificial_ops_,
         artificial_ops_ / dense_ops_ * 100);
  printf(" - est. max memory      ");
  printMemory(max_storage_);
}

void Symbolic::printShort() const {
  printf("\nStatistic of Factor L\n");
  printf("size            : %.2e\n", (double)n_);
  printf("nnz             : %.2e\n", nz_);
  printf("fill-in         : %.2f\n", fillin_);
  printf("operations      : %.2e\n", dense_ops_);
  printf("serial memory   : ");
  printMemory(max_storage_);
  printf("\n");
}