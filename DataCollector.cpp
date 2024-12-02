#include "DataCollector.h"

// instance of DataCollector
DataCollector* DataCollector::ptr_ = nullptr;

DataCollector::DataCollector() {
#ifdef DATA_COLLECTION
  times_.resize(kTimeSize);
  blas_calls_.resize(kTimeBlasEnd - kTimeBlasStart + 1);
#endif
}

DataCollector* DataCollector::get() {
  if (!ptr_) ptr_ = new DataCollector();
  return ptr_;
}

void DataCollector::destruct() {
  delete ptr_;
  ptr_ = nullptr;
}

void DataCollector::append() {
  // add an empty IterData object to the record
  iter_data_record_.push_back(IterData());
}

const IterData& DataCollector::iter(int i) const {
  // access the data record of a specific iteration
  return iter_data_record_[i];
}

IterData& DataCollector::last() {
  // access most recent record of data to write
  return iter_data_record_.back();
}
const IterData& DataCollector::last() const {
  // access most recent record of data to read
  return iter_data_record_.back();
}

void DataCollector::sumTime(TimeItems i, double t) {
#ifdef DATA_COLLECTION
  // Keep track of times and blas calls.
  std::lock_guard<std::mutex> lock(times_mutex_);
  times_[i] += t;
#ifdef BLAS_TIMING
  if (i >= kTimeBlasStart && i <= kTimeBlasEnd)
    ++blas_calls_[i - kTimeBlasStart];
#endif
#endif
}

void DataCollector::extremeEntries(double minD, double maxD, double minoffD,
                                   double maxoffD) {
#ifdef DATA_COLLECTION
  // Store max and min entries of D and L.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  last().minD = std::min(last().minD, minD);
  last().maxD = std::max(last().maxD, maxD);
  last().minL = std::min(last().minL, minoffD);
  last().maxL = std::max(last().maxL, maxoffD);
#endif
}

void DataCollector::sumRegPiv() {
#ifdef DATA_COLLECTION
  // Increase the number of dynamically regularized pivots.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++last().n_reg_piv;
#endif
}

void DataCollector::sumSwap() {
#ifdef DATA_COLLECTION
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++last().n_swap;
#endif
}

void DataCollector::sum2x2() {
#ifdef DATA_COLLECTION
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++last().n_2x2;
#endif
}

void DataCollector::setMaxReg(double new_reg) {
#ifdef DATA_COLLECTION
  // Keep track of maximum regularization used.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  last().max_reg = std::max(last().max_reg, new_reg);
#endif
}

void DataCollector::setWorstRes(double res) {
#ifdef DATA_COLLECTION
  // Keep track of worst residual
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  last().worst_res = std::max(last().worst_res, res);
#endif
}

void DataCollector::printTimes() const {
#ifdef COARSE_TIMING
  printf("----------------------------------------------------\n");
  printf("Analyse time            \t%8.4f\n", times_[kTimeAnalyse]);

#ifdef FINE_TIMING
  printf("\tMetis:                  %8.4f (%4.1f%%)\n",
         times_[kTimeAnalyseMetis],
         times_[kTimeAnalyseMetis] / times_[kTimeAnalyse] * 100);
  printf("\tTree:                   %8.4f (%4.1f%%)\n",
         times_[kTimeAnalyseTree],
         times_[kTimeAnalyseTree] / times_[kTimeAnalyse] * 100);
  printf("\tCounts:                 %8.4f (%4.1f%%)\n",
         times_[kTimeAnalyseCount],
         times_[kTimeAnalyseCount] / times_[kTimeAnalyse] * 100);
  printf("\tSupernodes:             %8.4f (%4.1f%%)\n", times_[kTimeAnalyseSn],
         times_[kTimeAnalyseSn] / times_[kTimeAnalyse] * 100);
  printf("\tReorder:                %8.4f (%4.1f%%)\n",
         times_[kTimeAnalyseReorder],
         times_[kTimeAnalyseReorder] / times_[kTimeAnalyse] * 100);
  printf("\tSn sparsity pattern:    %8.4f (%4.1f%%)\n",
         times_[kTimeAnalysePattern],
         times_[kTimeAnalysePattern] / times_[kTimeAnalyse] * 100);
  printf("\tRelative indices:       %8.4f (%4.1f%%)\n",
         times_[kTimeAnalyseRelInd],
         times_[kTimeAnalyseRelInd] / times_[kTimeAnalyse] * 100);
#endif

  printf("----------------------------------------------------\n");
  printf("Factorise time          \t%8.4f\n", times_[kTimeFactorise]);

#ifdef FINE_TIMING
  printf("\tPrepare:                %8.4f (%4.1f%%)\n",
         times_[kTimeFactorisePrepare],
         times_[kTimeFactorisePrepare] / times_[kTimeFactorise] * 100);
  printf("\tAssembly original:      %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseAssembleOriginal],
         times_[kTimeFactoriseAssembleOriginal] / times_[kTimeFactorise] * 100);
  printf("\tAssemble children in F: %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseAssembleChildrenFrontal],
         times_[kTimeFactoriseAssembleChildrenFrontal] /
             times_[kTimeFactorise] * 100);
  printf("\tAssemble children in C: %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseAssembleChildrenClique],
         times_[kTimeFactoriseAssembleChildrenClique] / times_[kTimeFactorise] *
             100);
  printf("\tDense factorisation:    %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseDenseFact],
         times_[kTimeFactoriseDenseFact] / times_[kTimeFactorise] * 100);
  printf("\t\tmain:           %8.4f\n", times_[kTimeDenseFact_main]);
  printf("\t\tSchur:          %8.4f\n", times_[kTimeDenseFact_schur]);
  printf("\t\tkernel:         %8.4f\n", times_[kTimeDenseFact_kernel]);
  printf("\t\tconvert:        %8.4f\n", times_[kTimeDenseFact_convert]);
  printf("\t\tpivoting:       %8.4f\n", times_[kTimeDenseFact_pivoting]);
  printf("\tTerminate:              %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseTerminate],
         times_[kTimeFactoriseTerminate] / times_[kTimeFactorise] * 100);
#endif

  printf("----------------------------------------------------\n");
  printf("Solve time              \t%8.4f\n", times_[kTimeSolve]);
  printf("----------------------------------------------------\n");

#ifdef BLAS_TIMING

  double total_blas_time =
      times_[kTimeBlas_copy] + times_[kTimeBlas_axpy] + times_[kTimeBlas_scal] +
      times_[kTimeBlas_swap] + times_[kTimeBlas_gemv] + times_[kTimeBlas_trsv] +
      times_[kTimeBlas_tpsv] + times_[kTimeBlas_ger] + times_[kTimeBlas_trsm] +
      times_[kTimeBlas_syrk] + times_[kTimeBlas_gemm];

  printf("BLAS time               \t%8.4f\n", total_blas_time);
  printf("\tcopy:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_copy], times_[kTimeBlas_copy] / total_blas_time * 100,
         blas_calls_[kTimeBlas_copy - kTimeBlasStart]);
  printf("\taxpy:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_axpy], times_[kTimeBlas_axpy] / total_blas_time * 100,
         blas_calls_[kTimeBlas_axpy - kTimeBlasStart]);
  printf("\tscal:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_scal], times_[kTimeBlas_scal] / total_blas_time * 100,
         blas_calls_[kTimeBlas_scal - kTimeBlasStart]);
  printf("\tswap:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_swap], times_[kTimeBlas_swap] / total_blas_time * 100,
         blas_calls_[kTimeBlas_swap - kTimeBlasStart]);
  printf("\tgemv:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_gemv], times_[kTimeBlas_gemv] / total_blas_time * 100,
         blas_calls_[kTimeBlas_gemv - kTimeBlasStart]);
  printf("\ttrsv:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_trsv], times_[kTimeBlas_trsv] / total_blas_time * 100,
         blas_calls_[kTimeBlas_trsv - kTimeBlasStart]);
  printf("\ttpsv:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_tpsv], times_[kTimeBlas_tpsv] / total_blas_time * 100,
         blas_calls_[kTimeBlas_tpsv - kTimeBlasStart]);
  printf("\tger:            \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_ger], times_[kTimeBlas_ger] / total_blas_time * 100,
         blas_calls_[kTimeBlas_ger - kTimeBlasStart]);
  printf("\ttrsm:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_trsm], times_[kTimeBlas_trsm] / total_blas_time * 100,
         blas_calls_[kTimeBlas_trsm - kTimeBlasStart]);
  printf("\tsyrk:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_syrk], times_[kTimeBlas_syrk] / total_blas_time * 100,
         blas_calls_[kTimeBlas_syrk - kTimeBlasStart]);
  printf("\tgemm:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_gemm], times_[kTimeBlas_gemm] / total_blas_time * 100,
         blas_calls_[kTimeBlas_gemm - kTimeBlasStart]);
  printf("----------------------------------------------------\n");
#endif
#endif
}
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
void DataCollector::printSymbolic(bool verbose) const {
  printf("\nStatistic of Factor L\n");
  printf("size            : %.2e\n", (double)n_);
  printf("nnz             : %.2e\n", nz_);
  printf("fill-in         : %.2f\n", fillin_);
  printf("serial memory   : ");
  printMemory(serial_storage_);

  printf("dense  ops      : %.1e\n", dense_ops_);
  printf("sparse ops      : %.1e\n", sparse_ops_);
  printf("critical ops    : %.1e\n", critical_ops_);
  printf("max tree speedup: %.2f\n", dense_ops_ / critical_ops_);

  if (verbose) {
    printf("artificial nz   : %.1e (%.1f%%)\n", (double)artificial_nz_,
           (double)artificial_nz_ / nz_ * 100);
    printf("artificial ops  : %.1e (%.1f%%)\n", artificial_ops_,
           artificial_ops_ / dense_ops_ * 100);
    printf("largest front   : %5d\n", largest_front_);
    printf("largest sn      : %5d\n", largest_sn_);
    printf("supernodes      : %5d\n", sn_);
    printf("sn size <=   1  : %5d\n", sn_size_1_);
    printf("sn size <=  10  : %5d\n", sn_size_10_);
    printf("sn size <= 100  : %5d\n", sn_size_100_);
  }
  printf("\n");
}

double DataCollector::minD() const { return last().minD; }
double DataCollector::maxD() const { return last().maxD; }
double DataCollector::minL() const { return last().minL; }
double DataCollector::maxL() const { return last().maxL; }
double DataCollector::maxReg() const { return last().max_reg; }
double DataCollector::worstRes() const { return last().worst_res; }
int DataCollector::nRegPiv() const { return last().n_reg_piv; }
int DataCollector::nSwaps() const { return last().n_swap; }
int DataCollector::n2x2() const { return last().n_2x2; }