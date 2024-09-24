#include "Symbolic.h"

#include <iostream>

FactType Symbolic::factType() const { return fact_type_; }
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

const std::vector<int>& Symbolic::ptr() const { return ptr_; }
const std::vector<int>& Symbolic::perm() const { return perm_; }
const std::vector<int>& Symbolic::iperm() const { return iperm_; }
const std::vector<int>& Symbolic::snParent() const { return sn_parent_; }
const std::vector<int>& Symbolic::snStart() const { return sn_start_; }
const std::vector<int>& Symbolic::pivotSign() const { return pivot_sign_; }

void Symbolic::setFact(FactType i) const { fact_type_ = i; }
void Symbolic::setFormat(FormatType i) const { format_type_ = i; }
double& Symbolic::times(TimeItems i) const { return times_record_[i]; }
std::vector<double>& Symbolic::dynamicReg() const { return dynamic_reg_; }
std::vector<double>& Symbolic::times() const { return times_record_; }

void Symbolic::print() const {
  printf("Symbolic factorisation:\n");
  printf(" - type                 %s\n",
         fact_type_ == FactType::Chol ? "Cholesky" : "LDLt");
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

void Symbolic::printTimes() const {
#ifdef COARSE_TIMING
  printf("----------------------------------------------------\n");
  printf("Analyse time            \t%8.4f\n", times_record_[kTimeAnalyse]);

#ifdef FINE_TIMING
  printf("\tMetis:                  %8.4f (%4.1f%%)\n",
         times_record_[kTimeAnalyseMetis],
         times_record_[kTimeAnalyseMetis] / times_record_[kTimeAnalyse] * 100);
  printf("\tTree:                   %8.4f (%4.1f%%)\n",
         times_record_[kTimeAnalyseTree],
         times_record_[kTimeAnalyseTree] / times_record_[kTimeAnalyse] * 100);
  printf("\tCounts:                 %8.4f (%4.1f%%)\n",
         times_record_[kTimeAnalyseCount],
         times_record_[kTimeAnalyseCount] / times_record_[kTimeAnalyse] * 100);
  printf("\tSupernodes:             %8.4f (%4.1f%%)\n",
         times_record_[kTimeAnalyseSn],
         times_record_[kTimeAnalyseSn] / times_record_[kTimeAnalyse] * 100);
  printf(
      "\tReorder:                %8.4f (%4.1f%%)\n",
      times_record_[kTimeAnalyseReorder],
      times_record_[kTimeAnalyseReorder] / times_record_[kTimeAnalyse] * 100);
  printf(
      "\tSn sparsity pattern:    %8.4f (%4.1f%%)\n",
      times_record_[kTimeAnalysePattern],
      times_record_[kTimeAnalysePattern] / times_record_[kTimeAnalyse] * 100);
  printf("\tRelative indices:       %8.4f (%4.1f%%)\n",
         times_record_[kTimeAnalyseRelInd],
         times_record_[kTimeAnalyseRelInd] / times_record_[kTimeAnalyse] * 100);
  printf("\tLayer 0:                %8.4f (%4.1f%%)\n",
         times_record_[kTimeAnalyseLayer0],
         times_record_[kTimeAnalyseLayer0] / times_record_[kTimeAnalyse] * 100);
#endif

  printf("----------------------------------------------------\n");
  printf("Factorise time          \t%8.4f\n", times_record_[kTimeFactorise]);

#ifdef FINE_TIMING
  printf("\tPrepare:                %8.4f (%4.1f%%)\n",
         times_record_[kTimeFactorisePrepare],
         times_record_[kTimeFactorisePrepare] / times_record_[kTimeFactorise] *
             100);
  printf("\tAssembly original:      %8.4f (%4.1f%%)\n",
         times_record_[kTimeFactoriseAssembleOriginal],
         times_record_[kTimeFactoriseAssembleOriginal] /
             times_record_[kTimeFactorise] * 100);
  printf("\tAssembly into frontal:  %8.4f (%4.1f%%)\n",
         times_record_[kTimeFactoriseAssembleChildrenF],
         times_record_[kTimeFactoriseAssembleChildrenF] /
             times_record_[kTimeFactorise] * 100);
  printf("\tAssembly into clique:   %8.4f (%4.1f%%)\n",
         times_record_[kTimeFactoriseAssembleChildrenC],
         times_record_[kTimeFactoriseAssembleChildrenC] /
             times_record_[kTimeFactorise] * 100);
  printf("\tDense factorisation:    %8.4f (%4.1f%%)\n",
         times_record_[kTimeFactoriseDenseFact],
         times_record_[kTimeFactoriseDenseFact] /
             times_record_[kTimeFactorise] * 100);
#endif

#ifdef FINEST_TIMING
  printf("\n");
  printf("\t\ttrsm:           %8.4f (%4.1f%%)\n",
         times_record_[kTimeDenseFact_trsm],
         times_record_[kTimeDenseFact_trsm] /
             times_record_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tsyrk:           %8.4f (%4.1f%%)\n",
         times_record_[kTimeDenseFact_syrk],
         times_record_[kTimeDenseFact_syrk] /
             times_record_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tgemm:           %8.4f (%4.1f%%)\n",
         times_record_[kTimeDenseFact_gemm],
         times_record_[kTimeDenseFact_gemm] /
             times_record_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tfact:           %8.4f (%4.1f%%)\n",
         times_record_[kTimeDenseFact_fact],
         times_record_[kTimeDenseFact_fact] /
             times_record_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tcopy:           %8.4f (%4.1f%%)\n",
         times_record_[kTimeDenseFact_copy],
         times_record_[kTimeDenseFact_copy] /
             times_record_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tcopy sch:       %8.4f (%4.1f%%)\n",
         times_record_[kTimeDenseFact_copyschur],
         times_record_[kTimeDenseFact_copyschur] /
             times_record_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tscal:           %8.4f (%4.1f%%)\n",
         times_record_[kTimeDenseFact_scal],
         times_record_[kTimeDenseFact_scal] /
             times_record_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tconvert:        %8.4f (%4.1f%%)\n",
         times_record_[kTimeDenseFact_convert],
         times_record_[kTimeDenseFact_convert] /
             times_record_[kTimeFactoriseDenseFact] * 100);
#endif
  printf("----------------------------------------------------\n");
  printf("Solve time              \t%8.4f\n", times_record_[kTimeSolve]);
  printf("----------------------------------------------------\n");

#endif
}