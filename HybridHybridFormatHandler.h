#ifndef HYBRID_HYBRID_FORMAT_H
#define HYBRID_HYBRID_FORMAT_H

#include "FormatHandler.h"

class HybridHybridFormatHandler : public FormatHandler {
  void initFrontal() override;
  int sizeClique() override;
  void assembleFrontal(int i, int j, double val) override;
  void assembleFrontalMultiple(int num, const double* child,
                               int nc, int child_sn, int row, int col, int i,
                               int j) override;
  int denseFactorise(double reg_thresh, std::vector<double>& regularization,int& n_reg_piv,
                     std::vector<double>& times) override;
  void assembleClique(const double* child, int nc,
                      int child_sn) override;
  void extremeEntries(double& minD, double& maxD, double& minoffD,
                      double& maxoffD) override;
};

#endif