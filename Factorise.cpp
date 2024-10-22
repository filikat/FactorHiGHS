#include "Factorise.h"

#include <fstream>

#include "FullFormatHandler.h"
#include "HybridHybridFormatHandler.h"
#include "HybridPackedFormatHandler.h"

Factorise::Factorise(const Symbolic& S, DataCollector& DC,
                     const std::vector<int>& rowsA,
                     const std::vector<int>& ptrA,
                     const std::vector<double>& valA)
    : S_{S}, DC_{DC} {
  // Input the symmetric matrix to be factorised in CSC format and the symbolic
  // factorisation coming from Analyze.
  // Only the lower triangular part of the matrix is used.

  n_ = ptrA.size() - 1;

  if (n_ != S_.size()) {
    printf(
        "Matrix provided to Factorise has size incompatible with symbolic "
        "object.\n");
    return;
  }

  rowsA_ = rowsA;
  valA_ = valA;
  ptrA_ = ptrA;

  // Permute the matrix.
  // This also removes any entry not in the lower triangle.
  permute(S_.iperm());

  nzA_ = ptrA_.back();

  // Double transpose to sort columns
  std::vector<int> temp_ptr(n_ + 1);
  std::vector<int> temp_rows(nzA_);
  std::vector<double> temp_val(nzA_);
  transpose(ptrA_, rowsA_, valA_, temp_ptr, temp_rows, temp_val);
  transpose(temp_ptr, temp_rows, temp_val, ptrA_, rowsA_, valA_);

  // create linked lists of children in supernodal elimination tree
  childrenLinkedList(S_.snParent(), first_children_, next_children_);

  // scale the matrix
  // equilibrate();

  // compute largest diagonal entry in absolute value
  max_diag_ = 0.0;
  min_diag_ = std::numeric_limits<double>::max();
  for (int col = 0; col < n_; ++col) {
    double val = std::fabs(valA_[ptrA_[col]]);
    max_diag_ = std::max(max_diag_, val);
    min_diag_ = std::min(min_diag_, val);
  }
  // printf("diag in [%e,%e]\n", min_diag_, max_diag_);
}

void Factorise::permute(const std::vector<int>& iperm) {
  // Symmetric permutation of the lower triangular matrix A based on inverse
  // permutation iperm.
  // The resulting matrix is lower triangular, regardless of the input matrix.

  std::vector<int> work(n_, 0);

  // go through the columns to count the nonzeros
  for (int j = 0; j < n_; ++j) {
    // get new index of column
    const int col = iperm[j];

    // go through elements of column
    for (int el = ptrA_[j]; el < ptrA_[j + 1]; ++el) {
      const int i = rowsA_[el];

      // ignore potential entries in upper triangular part
      if (i < j) continue;

      // get new index of row
      const int row = iperm[i];

      // since only lower triangular part is used, col is smaller than row
      int actual_col = std::min(row, col);
      ++work[actual_col];
    }
  }

  std::vector<int> new_ptr(n_ + 1);

  // get column pointers by summing the count of nonzeros in each column.
  // copy column pointers into work
  counts2Ptr(new_ptr, work);

  std::vector<int> new_rows(new_ptr.back());
  std::vector<double> new_val(new_ptr.back());

  // go through the columns to assign row indices
  for (int j = 0; j < n_; ++j) {
    // get new index of column
    const int col = iperm[j];

    // go through elements of column
    for (int el = ptrA_[j]; el < ptrA_[j + 1]; ++el) {
      const int i = rowsA_[el];

      // ignore potential entries in upper triangular part
      if (i < j) continue;

      // get new index of row
      const int row = iperm[i];

      // since only lower triangular part is used, col is smaller than row
      const int actual_col = std::min(row, col);
      const int actual_row = std::max(row, col);

      int pos = work[actual_col]++;
      new_rows[pos] = actual_row;
      new_val[pos] = valA_[el];
    }
  }

  ptrA_ = std::move(new_ptr);
  rowsA_ = std::move(new_rows);
  valA_ = std::move(new_val);
}

std::unique_ptr<FormatHandler> getFormatHandler(const Symbolic& S, int sn) {
  std::unique_ptr<FormatHandler> ptr;
  switch (S.formatType()) {
    case FormatType::Full:
      ptr.reset(new FullFormatHandler(S, sn));
      break;
    case FormatType::HybridPacked:
      ptr.reset(new HybridPackedFormatHandler(S, sn));
      break;
    case FormatType::HybridHybrid:
      ptr.reset(new HybridHybridFormatHandler(S, sn));
      break;
  }
  return ptr;
}

int Factorise::processSupernode(int sn) {
  // Assemble frontal matrix for supernode sn, perform partial factorisation and
  // store the result.
  Clock clock;

#ifdef FINE_TIMING
  clock.start();
#endif
  // ===================================================
  // Supernode information
  // ===================================================
  // first and last+1 column of the supernodes
  const int sn_begin = S_.snStart(sn);
  const int sn_end = S_.snStart(sn + 1);
  const int sn_size = sn_end - sn_begin;

  // initialize the format handler
  // this also allocates space for the frontal matrix and schur complement
  std::unique_ptr<FormatHandler> FH = getFormatHandler(S_, sn);

#ifdef FINE_TIMING
  DC_.times(kTimeFactorisePrepare) += clock.stop();
#endif

#ifdef FINE_TIMING
  clock.start();
#endif
  // ===================================================
  // Assemble original matrix A into frontal
  // ===================================================
  // j is relative column index in the frontal matrix
  for (int j = 0; j < sn_size; ++j) {
    // column index in the original matrix
    const int col = sn_begin + j;

    // go through the column
    for (int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
      // relative row index in the frontal matrix
      const int i = S_.relindCols(el);

      FH->assembleFrontal(i, j, valA_[el]);
    }
  }
#ifdef FINE_TIMING
  DC_.times(kTimeFactoriseAssembleOriginal) += clock.stop();
#endif

// ===================================================
// Assemble frontal matrices of children
// ===================================================
#ifdef FINE_TIMING
  clock.start();
#endif
  int child_sn = first_children_[sn];
  while (child_sn != -1) {
    // Schur contribution of the current child
    std::vector<double>& child_clique = schur_contribution_[child_sn];

    // determine size of clique of child
    const int child_begin = S_.snStart(child_sn);
    const int child_end = S_.snStart(child_sn + 1);

    // number of nodes in child sn
    const int child_size = child_end - child_begin;

    // size of clique of child sn
    const int nc = S_.ptr(child_sn + 1) - S_.ptr(child_sn) - child_size;

// ASSEMBLE INTO FRONTAL
#ifdef FINEST_TIMING
    Clock clock2;
    clock2.start();
#endif
    // go through the columns of the contribution of the child
    for (int col = 0; col < nc; ++col) {
      // relative index of column in the frontal matrix
      int j = S_.relindClique(child_sn, col);

      if (j < sn_size) {
        // assemble into frontal

        // go through the rows of the contribution of the child
        int row = col;
        while (row < nc) {
          // relative index of the entry in the matrix frontal
          const int i = S_.relindClique(child_sn, row);

          // how many entries to sum
          const int consecutive = S_.consecutiveSums(child_sn, row);

          FH->assembleFrontalMultiple(consecutive, child_clique, nc, child_sn,
                                      row, col, i, j);

          row += consecutive;
        }
      }
    }
#ifdef FINEST_TIMING
    DC_.times(kTimeFactoriseAssembleChildrenFrontal) += clock2.stop();
#endif

// ASSEMBLE INTO CLIQUE
#ifdef FINEST_TIMING
    clock2.start();
#endif
    FH->assembleClique(child_clique, nc, child_sn);
#ifdef FINEST_TIMING
    DC_.times(kTimeFactoriseAssembleChildrenClique) += clock2.stop();
#endif

    // Schur contribution of the child is no longer needed
    // Swap with temporary empty vector to deallocate memory
    std::vector<double> temp_empty;
    schur_contribution_[child_sn].swap(temp_empty);

    // move on to the next child
    child_sn = next_children_[child_sn];
  }
#ifdef FINE_TIMING
  DC_.times(kTimeFactoriseAssembleChildren) += clock.stop();
#endif

  // ===================================================
  // Partial factorisation
  // ===================================================
#ifdef FINE_TIMING
  clock.start();
#endif

  // threshold for regularization
  const double reg_thresh = max_diag_ * 1e-16 * 1e-10;

  int status = FH->denseFactorise(reg_thresh, DC_.n_reg_piv_, DC_.times());
  if (status) return status;

#ifdef FINE_TIMING
  DC_.times(kTimeFactoriseDenseFact) += clock.stop();
#endif

  // compute largest elements in factorization
  FH->extremeEntries(DC_);

  // terminate the format handler
  FH->terminate(sn_columns_[sn], schur_contribution_[sn], total_reg_);

  return kRetOk;
}

void Factorise::equilibrate() {
  // Scale the matrix, according to:
  // D. Ruiz, "A Scaling Algorithm to Equilibrate Both Rows and Columns Norms
  // in Matrices", RAL-TR-2001-034. It attempts to obtain a matrix with
  // /infty-norm of the rows and columns equal to one.

  const int maxiter = 10;

  // initialize equilibration factors
  colscale_.resize(n_, 1.0);

  // initialize vectors for max entry in cols
  std::vector<double> colmax(n_);

  // iterate
  for (int iter = 0; iter < maxiter; ++iter) {
    // compute \infty norm of cols
    colmax.assign(n_, 0.0);
    for (int col = 0; col < n_; ++col) {
      for (int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
        int row = rowsA_[el];
        double val = std::abs(valA_[el]);
        // Matrix is stored as lower triangle, so this element contributes
        // both to columns col and row.
        colmax[col] = std::max(colmax[col], val);
        if (col != row) colmax[row] = std::max(colmax[row], val);
      }
    }

    // check stopping criterion
    double max_deviation{};
    for (int i = 0; i < n_; ++i) {
      double val = std::abs(1.0 - colmax[i]);
      max_deviation = std::max(max_deviation, val);
    }
    if (max_deviation < .01) break;

    // compute scaling factors for columns
    for (int col = 0; col < n_; ++col) {
      // compute scaling
      colmax[col] = 1.0 / sqrt(colmax[col]);

      // round to power of 2
      int exp;
      std::frexp(colmax[col], &exp);
      colmax[col] = std::ldexp(1.0, exp);

      // update overall scaling
      colscale_[col] *= colmax[col];
    }

    // apply scaling to matrix
    for (int col = 0; col < n_; ++col) {
      for (int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
        int row = rowsA_[el];
        valA_[el] *= colmax[col];
        valA_[el] *= colmax[row];
      }
    }
  }

  /*double max_scale = 0.0;
  double min_scale = std::numeric_limits<double>::max();
  for (double d : colscale_) {
    max_scale = std::max(max_scale, d);
    min_scale = std::min(min_scale, d);
  }

  printf("Scaling in [%e,%e]\n", min_scale, max_scale);
  */
}

int Factorise::run(Numeric& num) {
  Clock clock;

#ifdef COARSE_TIMING
  clock.start();
#endif

  total_reg_.assign(n_, 0.0);

  // allocate space for list of generated elements and columns of L
  schur_contribution_.resize(S_.sn());
  sn_columns_.resize(S_.sn());

  DC_.resetExtremeEntries();

  int status{};
  for (int sn = 0; sn < S_.sn(); ++sn) {
    status = processSupernode(sn);
    if (status) break;
  }

  if (status) return status;

  // un-scale regularization
  if (colscale_.size() > 0) {
    for (int i = 0; i < n_; ++i) {
      total_reg_[i] /= (colscale_[i] * colscale_[i]);
    }
  }

  // compute maximum regularization
  DC_.max_reg_ = 0.0;
  for (int i = 0; i < total_reg_.size(); ++i) {
    DC_.max_reg_ = std::max(DC_.max_reg_, std::abs(total_reg_[i]));
  }

  // move factorisation to numerical object
  num.sn_columns_ = std::move(sn_columns_);
  num.colscale_ = std::move(colscale_);
  num.total_reg_ = std::move(total_reg_);

#ifdef COARSE_TIMING
  DC_.times(kTimeFactorise) += clock.stop();
#endif

  return kRetOk;
}