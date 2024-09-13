#ifndef HYBRID_PACKED_FORMAT_H
#define HYBRID_PACKED_FORMAT_H

#include "FormatHandler.h"

class HybridPackedFormatHandler : public FormatHandler {
  void initFrontal() override;
  void initClique() override;
  void assembleFrontal(int i, int j, double val) override;
  void assembleFrontalMultiple(int num, double* child, int nc, int child_sn,
                               int row, int col, int i, int j) override;
  int denseFactorise(std::vector<double>& times) override;
  void assembleClique(double* child, int nc, int child_sn) override;
};

#endif