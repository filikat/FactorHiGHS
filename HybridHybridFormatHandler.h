#ifndef HYBRID_HYBRID_FORMAT_H
#define HYBRID_HYBRID_FORMAT_H

#include "FormatHandler.h"

class HybridHybridFormatHandler : public FormatHandler {
  void initFrontal() override;
  void initClique() override;
  void assembleFrontal(int i, int j, double val) override;
  void assembleFrontalMultiple(int num, double* child, int nc, int child_sn,
                               int row, int col, int i, int j) override;
  int denseFactorise(double reg_thresh,std::vector<double>& regularization,std::vector<double>& times) override;
  void assembleClique(double* child, int nc, int child_sn) override;
};

#endif