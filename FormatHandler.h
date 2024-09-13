#ifndef FORMAT_HANDLER_H
#define FORMAT_HANDLER_H

#include <vector>

#include "Blas_declaration.h"
#include "DenseFact_declaration.h"
#include "Symbolic.h"

const int i_one = 1;
const double d_one = 1.0;

class FormatHandler {
 protected:
  std::vector<double>* frontal_{};
  double** clique_{};
  std::vector<std::vector<int>>* clique_block_start_{};
  const Symbolic* S_;

  int sn_{};
  int ldf_{};
  int ldc_{};
  int nb_{};
  int sn_size_{};

 public:
  void attach(std::vector<double>* frontal, double** clique,
              std::vector<std::vector<int>>* clique_block_start,
              const Symbolic* S, int sn);
  void detach();

  virtual void initFrontal() = 0;
  virtual void initClique() = 0;
  virtual void assembleFrontal(int i, int j, double val) = 0;
  virtual void assembleFrontalMultiple(int num, double* child, int nc,
                                       int child_sn, int row, int col, int i,
                                       int j) = 0;
  virtual int denseFactorise(std::vector<double>& times) = 0;
  virtual void assembleClique(double* child, int nc, int child_sn) = 0;
};

#endif