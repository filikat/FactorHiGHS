#include "Symbolic.h"

#include <iostream>

Symbolic::Symbolic(FactType fact_type, FormatType format_type, int n_threads)
    : fact_type_{fact_type}, format_type_{format_type}, n_threads_{n_threads} {}

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
int Symbolic::stackSize() const { return max_stack_entries_; }
int Symbolic::maxCliqueSize() const { return max_clique_entries_; }

const std::vector<int>& Symbolic::ptr() const { return ptr_; }
const std::vector<int>& Symbolic::iperm() const { return iperm_; }
const std::vector<int>& Symbolic::snParent() const { return sn_parent_; }
const std::vector<int>& Symbolic::snStart() const { return sn_start_; }
const std::vector<int>& Symbolic::pivotSign() const { return pivot_sign_; }

double& Symbolic::times(TimeItems i) const { return times_record_[i]; }
std::vector<double>& Symbolic::times() const { return times_record_; }

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
  printf(" - type                 %s\n",
         fact_type_ == FactType::Chol ? "Cholesky" : "LDLt");
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
  printf("max stack size  : %.2e\n", (double)max_stack_entries_);
  printf("\nRunning on %d thread%s\n", n_threads_, n_threads_ > 1 ? "s" : "");
  if (n_threads_ > 1) {
    printf("\nTree parallelization:\n");
    printf("thr  subtrees  total ops  ratio\n");

    const double max_load =
        *std::max_element(ops_per_thread_.begin(), ops_per_thread_.end());
    double total_ops = assembly_ops_ * 100 + dense_ops_;
    double ops_left = total_ops;

    for (int i = 0; i < n_threads_; ++i) {
      printf("%2d %7d %12.1e %7.2f\n", i, (int)subtrees_per_thread_[i].size(),
             ops_per_thread_[i], ops_per_thread_[i] / max_load);
      ops_left -= ops_per_thread_[i];
    }

    printf("Supernodes left  : %d\n", sn_above_layer0_);
    printf("Total ops left   : %.1e\n", ops_left);
    printf("Expected speedup : %.2f\n\n", total_ops / (ops_left + max_load));
  }
  printf("\n");
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
  printf("\t\taxpy:           %8.4f (%4.1f%%)\n",
         times_record_[kTimeDenseFact_axpy],
         times_record_[kTimeDenseFact_axpy] /
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