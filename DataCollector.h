#ifndef DATA_COLLECTOR_H
#define DATA_COLLECTOR_H

#include <atomic>
#include <mutex>
#include <vector>

#include "Timing.h"

struct IterData {
  // data of a given ipm iteration

  // factorization data
  double minD = std::numeric_limits<double>::max();
  double maxD = 0.0;
  double minL = std::numeric_limits<double>::max();
  double maxL = 0.0;
  double max_reg = 0.0;
  double worst_res = 0.0;
  int n_reg_piv = 0;
  int n_swap = 0;
  int n_2x2 = 0;
  int n_wrong_sign = 0;
  double max_wrong_sign = 0.0;

  // ipm data
  double p_obj;
  double d_obj;
  double p_inf;
  double d_inf;
  double mu;
  double pd_gap;
  double p_alpha;
  double d_alpha;
  double min_prod;
  double max_prod;
  double min_theta;
  double max_theta;
  double sigma;
  int correctors;
  int num_small_prod;
  int num_large_prod;
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

 public:
  // Access to the object
  static DataCollector* get();
  static void destruct();
  static IterData& back();

  // Manage record of data of iterations
  void append();

  // Functions with lock, they can be accessed simultaneously
  void sumTime(TimeItems i, double t);
  void countRegPiv();
  void countSwap();
  void count2x2();
  void setWrongSign(double p);
  void setMaxReg(double new_reg);
  void setWorstRes(double res);
  void setExtremeEntries(double minD, double maxD, double minoffD,
                         double maxoffD);

  // Const functions
  void printTimes() const;
  void printSymbolic(bool verbose = false) const;
  void printIter() const;
};

#endif