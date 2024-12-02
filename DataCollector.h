#ifndef DATA_COLLECTOR_H
#define DATA_COLLECTOR_H

#include <atomic>
#include <mutex>
#include <vector>

#include "Timing.h"

struct IterData {
  // data of a given ipm iteration
  double minD = std::numeric_limits<double>::max();
  double maxD = 0.0;
  double minL = std::numeric_limits<double>::max();
  double maxL = 0.0;
  double max_reg = 0.0;
  double worst_res = 0.0;
  int n_reg_piv = 0;
  int n_swap = 0;
  int n_2x2 = 0;
};

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
  int serial_storage_{};
  int largest_front_{};
  int largest_sn_{};
  int sn_size_1_{};
  int sn_size_10_{};
  int sn_size_100_{};

  // record of data of ipm iterations
  std::vector<IterData> iter_data_record_{};

  // Mutexes for concurrent access
  std::mutex times_mutex_;
  std::mutex iter_data_mutex_;

  // Instance of DataCollector
  static DataCollector* ptr_;

  friend class Analyse;

  // Private ctor and dtor
  DataCollector();
  ~DataCollector() = default;

  IterData& last();
  const IterData& last() const;

 public:
  // Access to the object
  static DataCollector* get();
  static void destruct();

  // Manage record of data of iterations
  void append();
  const IterData& iter(int i) const;

  // Functions with lock, they can be accessed simultaneously
  void sumTime(TimeItems i, double t);
  void sumRegPiv();
  void sumSwap();
  void sum2x2();
  void setMaxReg(double new_reg);
  void setWorstRes(double res);
  void extremeEntries(double minD, double maxD, double minoffD, double maxoffD);

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