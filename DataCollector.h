#ifndef DATA_COLLECTOR_H
#define DATA_COLLECTOR_H

#include <atomic>
#include <mutex>
#include <vector>

#include "Timing.h"

class DataCollector {
  // Record of times and BLAS calls
  std::vector<double> times_{};
  std::vector<int> blas_calls_{};

  // Symbolic factorization statistics
  int n_{};
  double nz_{};
  int sn_{};
  double fillin_{};
  double dense_ops_{};
  double sparse_ops_{};
  double critical_ops_{};
  int artificial_nz_{};
  double artificial_ops_{};
  int largest_front_{};
  int largest_sn_{};
  int serial_storage_{};

  // Other statistics
  double minD_{};
  double maxD_{};
  double minL_{};
  double maxL_{};
  double max_reg_{};
  double worst_res_{};
  int n_reg_piv_{};

  // Mutexes for concurrent access
  std::mutex times_mutex_;
  std::mutex extreme_entries_mutex_;
  std::mutex n_reg_piv_mutex_;
  std::mutex max_reg_mutex_;
  std::mutex worst_res_mutex_;

  friend class Analyse;

 public:
  // Constructor
  DataCollector();

  // Functions with lock, they can be accessed simultaneously
  void sumTime(TimeItems i, double t);
  void sumRegPiv();
  void setMaxReg(double new_reg);
  void setWorstRes(double res);
  void extremeEntries(double minD, double maxD, double minoffD, double maxoffD);

  // Functions without lock, to be accessed serially
  void resetExtremeEntries();

  // Const functions
  void printTimes() const;
  void printSymbolic(bool verbose = false) const;
  double minD() const;
  double maxD() const;
  double minL() const;
  double maxL() const;
  double maxReg() const;
  double worstRes() const;
  int nRegPiv() const;
};

#endif