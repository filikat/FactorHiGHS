#include "DataCollector.h"

DataCollector::DataCollector() { times_.resize(kTimeSize); }

void DataCollector::sumTime(TimeItems i, double t) {
  // will need lock
  times_[i] += t;
}

void DataCollector::resetExtremeEntries() {
  minD_ = std::numeric_limits<double>::max();
  maxD_ = 0.0;
  minL_ = std::numeric_limits<double>::max();
  maxL_ = 0.0;
  max_reg_ = 0.0;
  worst_res_ = 0.0;
  n_reg_piv_ = 0;
}

void DataCollector::extremeEntries(double minD, double maxD, double minoffD,
                                   double maxoffD) {
  // will need lock
  minD_ = std::min(minD_, minD);
  maxD_ = std::max(maxD_, maxD);
  minL_ = std::min(minL_, minoffD);
  maxL_ = std::max(maxL_, maxoffD);
}

void DataCollector::sumRegPiv() {
  // will need lock
  ++n_reg_piv_;
}
void DataCollector::setMaxReg(double new_reg) {
  // will need lock
  max_reg_ = std::max(max_reg_, new_reg);
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
  printf("\tAssembly children:      %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseAssembleChildren],
         times_[kTimeFactoriseAssembleChildren] / times_[kTimeFactorise] * 100);
#ifdef FINEST_TIMING
  printf("\t\tinto frontal:   %8.4f\n",
         times_[kTimeFactoriseAssembleChildrenFrontal]);
  printf("\t\tinto clique:    %8.4f\n",
         times_[kTimeFactoriseAssembleChildrenClique]);
#endif
  printf("\tDense factorisation:    %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseDenseFact],
         times_[kTimeFactoriseDenseFact] / times_[kTimeFactorise] * 100);
#endif

#ifdef FINEST_TIMING
  printf("\n");
  printf("\t\tcopy:           %8.4f (%4.1f%%)\n", times_[kTimeDenseFact_copy],
         times_[kTimeDenseFact_copy] / times_[kTimeFactoriseDenseFact] * 100);
  printf("\t\taxpy:           %8.4f (%4.1f%%)\n", times_[kTimeDenseFact_axpy],
         times_[kTimeDenseFact_axpy] / times_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tscal:           %8.4f (%4.1f%%)\n", times_[kTimeDenseFact_scal],
         times_[kTimeDenseFact_scal] / times_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tgemv:           %8.4f (%4.1f%%)\n", times_[kTimeDenseFact_gemv],
         times_[kTimeDenseFact_gemv] / times_[kTimeFactoriseDenseFact] * 100);
  printf("\t\ttrsm:           %8.4f (%4.1f%%)\n", times_[kTimeDenseFact_trsm],
         times_[kTimeDenseFact_trsm] / times_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tsyrk:           %8.4f (%4.1f%%)\n", times_[kTimeDenseFact_syrk],
         times_[kTimeDenseFact_syrk] / times_[kTimeFactoriseDenseFact] * 100);
  printf("\t\tgemm:           %8.4f (%4.1f%%)\n", times_[kTimeDenseFact_gemm],
         times_[kTimeDenseFact_gemm] / times_[kTimeFactoriseDenseFact] * 100);

  printf("\t\tfact:           %8.4f\n", times_[kTimeDenseFact_fact]);
  printf("\t\tconvert:        %8.4f\n", times_[kTimeDenseFact_convert]);
#endif
  printf("----------------------------------------------------\n");
  printf("Solve time              \t%8.4f\n", times_[kTimeSolve]);
  printf("----------------------------------------------------\n");

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
  printf("operations      : %.2e\n", dense_ops_);
  printf("serial memory   : ");
  printMemory(serial_storage_);
  if (verbose) {
    printf("sparse ops      : %.2e\n", sparse_ops_);
    printf("artificial nz   : %.2e\n", (double)artificial_nz_);
    printf("artificial ops  : %.2e\n", artificial_ops_);
    printf("largest front   : %.2e\n", (double)largest_front_);
    printf("largest sn      : %.2e\n", (double)largest_sn_);
  }
  printf("\n");
}