#ifndef HYBRID_SOLVE_HANDLER_H
#define HYBRID_SOLVE_HANDLER_H

#include "Auxiliary.h"
#include "SolveHandler.h"

class HybridSolveHandler : public SolveHandler {
  void forwardSolve(std::vector<double>& x) const override;
  void backwardSolve(std::vector<double>& x) const override;
  void diagSolve(std::vector<double>& x) const override;

 public:
  HybridSolveHandler(const Symbolic& S, DataCollector& DC,
                     const std::vector<std::vector<double>>& sn_columns);
};

#endif