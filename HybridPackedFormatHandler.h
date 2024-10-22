#ifndef HYBRID_PACKED_FORMAT_H
#define HYBRID_PACKED_FORMAT_H

#include "FormatHandler.h"

class HybridPackedFormatHandler : public FormatHandler {
  void initFrontal() override;
  void initClique() override;
  void assembleFrontal(int i, int j, double val) override;
  void assembleFrontalMultiple(int num, const std::vector<double>& child,
                               int nc, int child_sn, int row, int col, int i,
                               int j) override;
  int denseFactorise(double reg_thresh, std::vector<double>& regularization,
                     int& n_reg_piv, std::vector<double>& times) override;
  void assembleClique(const std::vector<double>& child, int nc,
                      int child_sn) override;
  void extremeEntries(DataCollector& DC) override;

 public:
  HybridPackedFormatHandler(const Symbolic& S, int sn);
};

#endif