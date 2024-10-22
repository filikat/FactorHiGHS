#ifndef DATA_COLLECTOR_H
#define DATA_COLLECTOR_H

#include <vector>

#include "timing.h"

class DataCollector {
 public:
  // ==== Times ====
  std::vector<double> times_{};

  double& times(TimeItems i);
  std::vector<double>& times();

  void printTimes() const;

  // ==== Symbolic factorization statistics ====
  int n_{};
  double nz_{};
  double fillin_{};
  double dense_ops_{};
  double sparse_ops_{};
  int artificial_nz_{};
  double artificial_ops_{};
  int largest_front_{};
  int largest_sn_{};
  int serial_storage_{};
  void printSymbolic(bool verbose = false) const;

  // ==== Data ====
  double minD_, maxD_, minL_, maxL_;
  double max_reg_;
  double worst_res_;
  int n_reg_piv_;

  void resetExtremeEntries();
  void extremeEntries(double minD, double maxD, double minoffD, double maxoffD);
};

#endif