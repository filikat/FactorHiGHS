#ifndef FULL_FORMAT_H
#define FULL_FORMAT_H

#include "FormatHandler.h"

class FullFormatHandler : public FormatHandler {
  void initFrontal() override;
  void initClique() override;
  void assembleFrontal(int i, int j, double val) override;
  void assembleFrontalMultiple(int num, const std::vector<double>& child,
                               int nc, int child_sn, int row, int col, int i,
                               int j) override;
  int denseFactorise(double reg_thresh, std::vector<double>& regularization,
                     std::vector<double>& times) override;
  void assembleClique(const std::vector<double>& child, int nc,
                      int child_sn) override;
  void extremeEntries(double& minD, double& maxD, double& minoffD,
                      double& maxoffD) override;
};

#endif