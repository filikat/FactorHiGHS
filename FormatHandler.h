#ifndef FORMAT_HANDLER_H
#define FORMAT_HANDLER_H

#include <vector>

#include "Auxiliary.h"
#include "Blas_declaration.h"
#include "DataCollector.h"
#include "DenseFact_declaration.h"
#include "Symbolic.h"

// Interface class to handle different formats of dense matrices during the
// factorise phase.
// Any implementation of a specific format needs to define:
// - initFrontal: to initialize the frontal matrix with the correct number of
//                elements; the entries should be set to zero.
// - initClique: to initialize the clique matrix with the correct number of
//                elements; the entries should be left uninitialized.
// - assembleFrontal: to set a specific entry of frontal (used to assemble the
//                original matrix)
// - assembleFrontalMultiple: to sum a given number of consecutive entries into
//                frontal (used to assemble the child supernodes)
// - denseFactorise: to perform the dense partial factorization of frontal,
//                storing the Schur complement in clique.
// - assembleClique: to sum the contributions of a given child supernode into
//                clique.
//

// Constants for BLAS calls
const int i_one = 1;
const double d_one = 1.0;

class FormatHandler {
 protected:
  // data shared by all supernodes
  std::vector<std::vector<int>> clique_block_start_{};
  const Symbolic* S_;

  // data of a given supernode
  std::vector<double> frontal_{};
  std::vector<double> clique_{};
  std::vector<int> clique_block_start_sn_{};

  // which supernode is being processed
  int sn_{};

  // size of the front
  int ldf_{};

  // size of the clique
  int ldc_{};

  // block size
  int nb_{};

  // size of the supernode
  int sn_size_{};

 public:
  // initialize the whole object
  void init(const Symbolic* S);

  // initialize the data of a specific supernode
  void attach(int sn);

  // reset the FormatHandler
  void detach(std::vector<double>& frontal, std::vector<double>& clique);

  // =================================================================
  // Pure virtual functions.
  // These need to be defined by any derived class.
  // =================================================================
  virtual void initFrontal() = 0;
  virtual void initClique() = 0;
  virtual void assembleFrontal(int i, int j, double val) = 0;
  virtual void assembleFrontalMultiple(int num,
                                       const std::vector<double>& child, int nc,
                                       int child_sn, int row, int col, int i,
                                       int j) = 0;
  virtual int denseFactorise(double reg_thresh,
                             std::vector<double>& regularization,
                             int& n_reg_piv, std::vector<double>& times) = 0;
  virtual void assembleClique(const std::vector<double>& child, int nc,
                              int child_sn) = 0;
  // =================================================================

  // =================================================================
  // Virtual functions.
  // These may be overridden by derived classes, if needed.
  // =================================================================
  virtual void extremeEntries(DataCollector& DC) {}
  // =================================================================
};

#endif