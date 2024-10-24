#include "Factorise.h"

#include <fstream>

#include "FullFormatHandler.h"
#include "HybridHybridFormatHandler.h"
#include "HybridPackedFormatHandler.h"

Factorise::Factorise(const Symbolic& S, CliqueStack& stack,
                     const std::vector<int>& rowsA,
                     const std::vector<int>& ptrA,
                     const std::vector<double>& valA)
    : S_{S}, clique_stack_{stack} {
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

  // scale();

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

  // Allocate space for frontal matrix:
  // The front size is ldf and the supernode size is sn_size.
  // The frontal matrix is stored as two dense matrices:
  // frontal is ldf x sn_size and stores the first sn_size columns that will
  // undergo Cholesky elimination.
  // clique is ldc x ldc and stores the remaining (ldf - sn_size) columns
  // (without the top part), that do not undergo Cholesky elimination.

  // attach to the format handler
  FH_->attach(&sn_columns_[sn], sn);

  // initialize frontal
  FH_->initFrontal();

  // initialize clique
  int clique_size = FH_->sizeClique();
  double* clique = clique_stack_.setup(clique_size);
  FH_->attachClique(clique);

#ifdef FINE_TIMING
  S_.times(kTimeFactorisePrepare) += clock.stop();
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

      FH_->assembleFrontal(i, j, valA_[el]);
    }
  }
#ifdef FINE_TIMING
  S_.times(kTimeFactoriseAssembleOriginal) += clock.stop();
#endif

  // ===================================================
  // Assemble frontal matrices of children
  // ===================================================
  // child are assembled in reverse order (because they are stored in a stack).
  // we still use the linked list to iterate through them, to get the right
  // number of children, but the child id needs to be read from the stack.

  int child = first_children_[sn];
  while (child != -1) {
    // Schur contribution of the current child
    int child_sn{};
    const double* child_clique = clique_stack_.getChild(child_sn);

    // determine size of clique of child
    const int child_begin = S_.snStart(child_sn);
    const int child_end = S_.snStart(child_sn + 1);

    // number of nodes in child sn
    const int child_size = child_end - child_begin;

    // size of clique of child sn
    const int nc = S_.ptr(child_sn + 1) - S_.ptr(child_sn) - child_size;

// ASSEMBLE INTO FRONTAL
#ifdef FINE_TIMING
    clock.start();
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

          FH_->assembleFrontalMultiple(consecutive, child_clique, nc, child_sn,
                                       row, col, i, j);

          row += consecutive;
        }
      }
    }
#ifdef FINE_TIMING
    S_.times(kTimeFactoriseAssembleChildrenF) += clock.stop();
#endif

// ASSEMBLE INTO CLIQUE
#ifdef FINE_TIMING
    clock.start();
#endif
    FH_->assembleClique(child_clique, nc, child_sn);
#ifdef FINE_TIMING
    S_.times(kTimeFactoriseAssembleChildrenC) += clock.stop();
#endif

    // child is no longer needed, remove it from stack
    clique_stack_.popChild();

    // move on to the next child
    child = next_children_[child];
  }

  // ===================================================
  // Partial factorisation
  // ===================================================
#ifdef FINE_TIMING
  clock.start();
#endif

  // threshold for regularization
  const double reg_thresh = max_diag_ * 1e-16 * 1e-10;

  int status =
      FH_->denseFactorise(reg_thresh, total_reg_, n_reg_piv_, S_.times());
  if (status) return status;

#ifdef FINE_TIMING
  S_.times(kTimeFactoriseDenseFact) += clock.stop();
#endif

  // compute largest elements in factorization
  double minD{}, maxD{}, minoffD{}, maxoffD{};
  FH_->extremeEntries(minD, maxD, minoffD, maxoffD);
  minD_ = std::min(minD_, minD);
  maxD_ = std::max(maxD_, maxD);
  minoffD_ = std::min(minoffD_, minoffD);
  maxoffD_ = std::max(maxoffD_, maxoffD);

  // detach from the format handler
  FH_->detach();

  clique_stack_.pushWork(sn);

  return kRetOk;
}

bool Factorise::check() const {
  // Check that the numerical factorisation is correct, by using dense linear
  // algebra operations.
  // Return true if check is successful, or if matrix is too large.
  // To be used for debug.

  if (S_.factType() == FactType::LDLt ||
      S_.formatType() == FormatType::HybridPacked) {
    printf("\n==> Dense check not available\n");
    return true;
  }

  if (n_ > 5000) {
    printf("\n==> Matrix is too large for dense check\n\n");
    return true;
  }

  // assemble sparse matrix into dense matrix
  std::vector<double> M(n_ * n_);
  for (int col = 0; col < n_; ++col) {
    for (int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
      int row = rowsA_[el];

      // insert element in position (row,col)
      M[row + col * n_] = valA_[el];
    }
  }

  // use Lapack to factorise the dense matrix
  char uplo = 'L';
  int N = n_;
  int info;
  dpotrf_(&uplo, &N, M.data(), &N, &info);
  if (info != 0) {
    printf("\n==> dpotrf failed\n\n");
    return false;
  }

  // assemble sparse factor into dense factor
  std::vector<double> L(n_ * n_);
  for (int sn = 0; sn < S_.sn(); ++sn) {
    for (int col = S_.snStart(sn); col < S_.snStart(sn + 1); ++col) {
      // indices to access corresponding entry in the supernode
      int col_sn = col - S_.snStart(sn);
      int ldsn = S_.ptr(sn + 1) - S_.ptr(sn);

      for (int el = S_.ptr(sn); el < S_.ptr(sn + 1); ++el) {
        int row = S_.rows(el);

        // indices to access corresponding entry in the supernode
        int row_sn = el - S_.ptr(sn);

        // skip upper triangle of supernodes
        if (row < col) continue;

        L[row + col * n_] = sn_columns_[sn][row_sn + col_sn * ldsn];
      }
    }
  }

  // Check that sparse factorisation agrees with dense one.
  // This is done by computing the Frobenius norm of the difference between
  // the dense and sparse factors, divided by the Frobenius norm of the dense
  // factor.

  double frobenius_dense{};
  double frobenius_diff{};

  for (int col = 0; col < n_; ++col) {
    for (int row = 0; row < n_; ++row) {
      double val_sparse = L[row + n_ * col];
      double val_dense = M[row + col * n_];
      double diff = val_sparse - val_dense;

      frobenius_dense += val_dense * val_dense;
      frobenius_diff += diff * diff;
    }
  }

  frobenius_dense = sqrt(frobenius_dense);
  frobenius_diff = sqrt(frobenius_diff);
  double check_error = frobenius_diff / frobenius_dense;

  printf("\nFactorise Frobenius error %e\n", check_error);
  if (check_error < 1e-12) {
    printf("\n==> Factorise check successful\n\n");
    return true;
  } else {
    printf("\n==> Factorise check failed\n\n");
    return false;
  }
}

void Factorise::equilibrate() {
  // Scale the matrix, according to:
  // D. Ruiz, "A Scaling Algorithm to Equilibrate Both Rows and Columns Norms in
  // Matrices", RAL-TR-2001-034.
  // It attempts to obtain a matrix with /infty-norm of the rows and columns
  // equal to one.

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
        // Matrix is stored as lower triangle, so this element contributes both
        // to columns col and row.
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

void Factorise::scale() {
  colexp_.resize(n_);
  CurtisReidScalingSym(ptrA_, rowsA_, valA_, colexp_);

  for (int col = 0; col < n_; ++col) {
    for (int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
      int row = rowsA_[el];
      valA_[el] = std::ldexp(valA_[el], colexp_[col]);
      valA_[el] = std::ldexp(valA_[el], colexp_[row]);
    }
  }

  int maxexp = 0;
  int minexp = INT_MAX;
  for (int i = 0; i < n_; ++i) {
    maxexp = std::max(maxexp, colexp_[i]);
    minexp = std::min(minexp, colexp_[i]);
  }
  printf("Scaling: [%.2e,%.2e]\n", std::ldexp(1.0, minexp),
         std::ldexp(1.0, maxexp));
}

int Factorise::run(Numeric& num) {
  Clock clock;

#ifdef COARSE_TIMING
  clock.start();
#endif

  total_reg_.assign(n_, 0.0);

  // Handle multiple formats
  FullFormatHandler full_FH;
  HybridPackedFormatHandler hybrid_packed_FH;
  HybridHybridFormatHandler hybrid_hybrid_FH;

  switch (S_.formatType()) {
    case FormatType::Full:
      FH_ = &full_FH;
      break;
    case FormatType::HybridPacked:
      FH_ = &hybrid_packed_FH;
      break;
    case FormatType::HybridHybrid:
      FH_ = &hybrid_hybrid_FH;
      break;
  }

  // initialize format handler
  FH_->init(&S_);

  // allocate space for columns of L
  sn_columns_.resize(S_.sn());

  // initialize stack of cliques
  clique_stack_.init(S_.stackSize(), S_.maxCliqueSize());

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
  } else if (colexp_.size() > 0) {
    for (int i = 0; i < n_; ++i) {
      total_reg_[i] = std::ldexp(total_reg_[i], -2 * colexp_[i]);
    }
  }

  // move factorisation to numerical object
  num.sn_columns_ = std::move(sn_columns_);
  num.S_ = &S_;
  num.colscale_ = std::move(colscale_);
  num.colexp_ = std::move(colexp_);
  num.total_reg_ = std::move(total_reg_);
  num.n_reg_piv_ = n_reg_piv_;

#ifdef COARSE_TIMING
  S_.times(kTimeFactorise) += clock.stop();
#endif

  return kRetOk;
}