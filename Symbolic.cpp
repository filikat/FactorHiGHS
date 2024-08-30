#include "Symbolic.h"

#include <iostream>

void Symbolic::print() const {
  printf("Symbolic factorisation:\n");
  printf(" - type                 %s\n",
         type_ == FactType::NormEq ? "Normal equations" : "Augmented system");
  printf(" - size                 %d\n", n_);
  printf(" - nonzero entries      %.2e\n", (double)nz_);
  printf(" - density              %.2f\n", ((double)nz_ / n_) / n_);
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

  if (max_storage_ > 0) {
    printf(" - est. max memory      ");
    if (max_storage_ < 1024) {
      printf("%.2f Bytes\n", max_storage_);
    } else if (max_storage_ < 1024 * 1024) {
      printf("%.2f KB\n", max_storage_ / 1024);
    } else if (max_storage_ < 1024 * 1024 * 1024) {
      printf("%.2f MB\n", max_storage_ / 1024 / 1024);
    } else {
      printf("%.2f GB\n", max_storage_ / 1024 / 1024 / 1024);
    }
  }
}

FactType Symbolic::type() const { return type_; }
PackType Symbolic::packFormat() const { return pack_format_; }
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

const std::vector<int>& Symbolic::ptr() const { return ptr_; }
const std::vector<int>& Symbolic::perm() const { return perm_; }
const std::vector<int>& Symbolic::iperm() const { return iperm_; }
const std::vector<int>& Symbolic::snParent() const { return sn_parent_; }
const std::vector<int>& Symbolic::snStart() const { return sn_start_; }
