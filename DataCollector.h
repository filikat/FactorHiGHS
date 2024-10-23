#ifndef DATA_COLLECTOR_H
#define DATA_COLLECTOR_H

#include <vector>

#include "Timing.h"

class DataCollector {
 public:
  DataCollector();

  // ==== Times ====
  std::vector<double> times_{};

  void sumTime(TimeItems i, double t);
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

  // ==== Other statistics ====
  double minD_, maxD_, minL_, maxL_;
  double max_reg_;
  double worst_res_;
  int n_reg_piv_;

  void sumRegPiv();
  void setMaxReg(double new_reg);

  void resetExtremeEntries();
  void extremeEntries(double minD, double maxD, double minoffD, double maxoffD);
};

#endif