#ifndef DATA_COLLECTOR_H
#define DATA_COLLECTOR_H

#include <atomic>
#include <mutex>
#include <vector>

#include "Timing.h"

// DataCollector is a singleton object.
// Only one copy of it can exist and it cannot be constructed or destructed
// explicitly. Any public member function should be accessed through
// DataCollector::get()-> ...

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
  std::atomic<int> n_reg_piv_{};
  std::atomic<int> n_swaps_{};
  std::atomic<int> n_2x2_{};

  // Mutexes for concurrent access
  std::mutex times_mutex_;
  std::mutex extreme_entries_mutex_;
  std::mutex max_reg_mutex_;
  std::mutex worst_res_mutex_;

  // Instance of DataCollector
  static DataCollector* ptr_;

  friend class Analyse;

  // Private ctor and dtor
  DataCollector();
  ~DataCollector() = default;

 public:
  // Access to the object
  static DataCollector* get();
  static void destruct();

  // Functions with lock, they can be accessed simultaneously
  void sumTime(TimeItems i, double t);
  void sumRegPiv();
  void sumSwap();
  void sum2x2();
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
  int nSwaps() const;
  int n2x2() const;
};

#endif