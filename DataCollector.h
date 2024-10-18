#ifndef DATA_COLLECTOR_H
#define DATA_COLLECTOR_H

#include <vector>

#include "timing.h"

class DataCollector {
  // Times
 private:
  std::vector<double> times_{};

 public:
  double& times(TimeItems i);
  std::vector<double>& times();

  void printTimes() const;

  // Data
 public:
  double minD_, maxD_, minL_, maxL_;
  int num_reg_;
  double max_reg_;
  double worst_res_;
  int n_reg_piv_;

  void resetExtremeEntries();
  void extremeEntries(double minD, double maxD, double minoffD, double maxoffD);
};

#endif