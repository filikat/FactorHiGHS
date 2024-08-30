#include "Factorise.h"

#include <fstream>

Factorise::Factorise(const Symbolic& S, const std::vector<int>& rowsA,
                     const std::vector<int>& ptrA,
                     const std::vector<double>& valA)
    : S_{S} {
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
  // This also removes any entry not in the upper triangle.
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

  // allocate space for list of generated elements and columns of L
  schur_contribution_.resize(S_.sn(), nullptr);
  sn_columns_.resize(S_.sn());
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

  clock.start();
  // ===================================================
  // Supernode information
  // ===================================================
  // first and last+1 column of the supernodes
  const int sn_begin = S_.snStart(sn);
  const int sn_end = S_.snStart(sn + 1);
  const int sn_size = sn_end - sn_begin;

  // leading dimension of the frontal matrix
  const int ldf = S_.ptr(sn + 1) - S_.ptr(sn);

  // leading dimension of the clique matrix
  const int ldc = ldf - sn_size;

  // Allocate space for frontal matrix:
  // The front size is ldf and the supernode size is sn_size.
  // The frontal matrix is stored as two dense matrices:
  // frontal is ldf x sn_size and stores the first sn_size columns that will
  // undergo Cholesky elimination.
  // clique is ldc x ldc and stores the remaining (ldf - sn_size) columns
  // (without the top part), that do not undergo Cholesky elimination.
  std::vector<double>& frontal = sn_columns_[sn];
  double*& clique = schur_contribution_[sn];

  // frontal is initialized to zero
  switch (S_.packFormat()) {
    case PackType::Full:
      frontal.resize(ldf * sn_size, 0.0);
      break;
    case PackType::Hybrid:
    case PackType::Hybrid2:
      frontal.resize(ldf * sn_size - sn_size * (sn_size - 1) / 2, 0.0);
      break;
  }

  // clique need not be initialized to zero, provided that the assembly is done
  // properly
  switch (S_.packFormat()) {
    case PackType::Full:
      if (ldc > 0) clique = new double[ldc * ldc];
      break;

    case PackType::Hybrid2:
    case PackType::Hybrid: {
      const int nb = S_.blockSize();
      const int n_blocks = (ldc - 1) / nb + 1;
      clique_block_start_[sn].resize(n_blocks + 1);
      int schur_size{};
      for (int j = 0; j < n_blocks; ++j) {
        clique_block_start_[sn][j] = schur_size;
        const int jb = std::min(nb, ldc - j * nb);
        schur_size += (ldc - j * nb) * jb;
      }
      clique_block_start_[sn].back() = schur_size;
      clique = new double[schur_size];
    } break;
  }

  time_prepare_ += clock.stop();

  clock.start();
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

      switch (S_.packFormat()) {
        case PackType::Full:
          frontal[i + j * ldf] = valA_[el];
          break;
        case PackType::Hybrid:
        case PackType::Hybrid2:
          frontal[i + j * ldf - j * (j + 1) / 2] = valA_[el];
          break;
      }
    }
  }
  time_assemble_original_ += clock.stop();

  // ===================================================
  // Assemble frontal matrices of children into frontal
  // ===================================================
  clock.start();
  int child_sn = first_children_[sn];
  while (child_sn != -1) {
    // Schur contribution of the current child
    double* child_clique = schur_contribution_[child_sn];
    if (!child_clique) {
      printf("Error with child supernode\n");
      return kRetGeneric;
    }

    // determine size of clique of child
    const int child_begin = S_.snStart(child_sn);
    const int child_end = S_.snStart(child_sn + 1);

    // number of nodes in child sn
    const int child_size = child_end - child_begin;

    // size of clique of child sn
    const int nc = S_.ptr(child_sn + 1) - S_.ptr(child_sn) - child_size;

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

          // use daxpy_ for summing consecutive entries
          const int i_one = 1;
          const double d_one = 1.0;
          switch (S_.packFormat()) {
            case PackType::Full:
              daxpy_(&consecutive, &d_one, &child_clique[row + nc * col],
                     &i_one, &frontal[i + ldf * j], &i_one);
              break;

            case PackType::Hybrid2: {
              const int nb = S_.blockSize();
              const int jblock = col / nb;
              const int jb = std::min(nb, nc - nb * jblock);
              const int row_ = row - jblock * nb;
              const int col_ = col - jblock * nb;
              const int start_block = clique_block_start_[child_sn][jblock];
              daxpy_(&consecutive, &d_one,
                     &child_clique[start_block + col_ + jb * row_], &jb,
                     &frontal[i + ldf * j - j * (j + 1) / 2], &i_one);
            } break;

            case PackType::Hybrid: {
              const int nb = S_.blockSize();
              const int jblock = col / nb;
              const int row_ = row - jblock * nb;
              const int col_ = col - jblock * nb;
              const int start_block = clique_block_start_[child_sn][jblock];
              const int ld = nc - nb * jblock;
              daxpy_(&consecutive, &d_one,
                     &child_clique[start_block + row_ + ld * col_], &i_one,
                     &frontal[i + ldf * j - j * (j + 1) / 2], &i_one);
            } break;
          }
          row += consecutive;
        }
      }

      // If j >= sn_size, we would assemble into clique.
      // This is delayed until after the partial factorisation, to avoid
      // having to initialize clique to zero.
    }

    // move on to the next child
    child_sn = next_children_[child_sn];
  }
  time_assemble_children_F_ += clock.stop();

  // ===================================================
  // Partial factorisation
  // ===================================================
  clock.start();
  switch (S_.packFormat()) {
    case PackType::Full:
      if (S_.type() == FactType::NormEq) {
        int status =
            dense_fact_pdbf(ldf, sn_size, S_.blockSize(), frontal.data(), ldf,
                           clique, ldc, times_dense_fact_.data());
        if (status) return status;

      } else {
        int status =
            dense_fact_pibf(ldf, sn_size, S_.blockSize(), frontal.data(), ldf,
                           clique, ldc, times_dense_fact_.data());
        if (status) return status;
      }
      break;

    case PackType::Hybrid2: {
      int status = dense_fact_l2h(frontal.data(), ldf, sn_size, S_.blockSize(),
                                 times_dense_fact_.data());
      if (status) return status;

      if (S_.type() == FactType::NormEq) {
        status = dense_fact_pdbh_2(ldf, sn_size, S_.blockSize(), frontal.data(),
                                  clique, times_dense_fact_.data());
        if (status) return status;
      } else {
        status = dense_fact_pibh_2(ldf, sn_size, S_.blockSize(), frontal.data(),
                                  clique, times_dense_fact_.data());
        if (status) return status;
      }
    } break;

    case PackType::Hybrid: {
      int status = dense_fact_l2h(frontal.data(), ldf, sn_size, S_.blockSize(),
                                 times_dense_fact_.data());
      if (status) return status;

      if (S_.type() == FactType::NormEq) {
        status = dense_fact_pdbh(ldf, sn_size, S_.blockSize(), frontal.data(),
                                clique, times_dense_fact_.data());
        if (status) return status;
      } else {
        status = dense_fact_pibh(ldf, sn_size, S_.blockSize(), frontal.data(),
                                clique, times_dense_fact_.data());
        if (status) return status;
      }
    } break;
  }

  time_factorise_ += clock.stop();

  // ===================================================
  // Assemble frontal matrices of children into clique
  // ===================================================
  clock.start();
  child_sn = first_children_[sn];
  while (child_sn != -1) {
    // Schur contribution of the current child
    double* child_clique = schur_contribution_[child_sn];
    if (!child_clique) {
      printf("Error with child supernode\n");
      return kRetGeneric;
    }

    // determine size of clique of child
    const int child_begin = S_.snStart(child_sn);
    const int child_end = S_.snStart(child_sn + 1);

    // number of nodes in child sn
    const int child_size = child_end - child_begin;

    // size of clique of child sn
    const int nc = S_.ptr(child_sn + 1) - S_.ptr(child_sn) - child_size;

    if (S_.packFormat() != PackType::Hybrid2) {
      //   if (true) {
      //   go through the columns of the contribution of the child
      for (int col = 0; col < nc; ++col) {
        // relative index of column in the frontal matrix
        int j = S_.relindClique(child_sn, col);

        if (j >= sn_size) {
          // assemble into clique

          // adjust relative index to access clique
          j -= sn_size;

          // go through the rows of the contribution of the child
          int row = col;
          while (row < nc) {
            // relative index of the entry in the matrix clique
            const int i = S_.relindClique(child_sn, row) - sn_size;

            // how many entries to sum
            const int consecutive = S_.consecutiveSums(child_sn, row);

            // use daxpy_ for summing consecutive entries
            const int i_one = 1;
            const double d_one = 1.0;
            switch (S_.packFormat()) {
              case PackType::Full:
                daxpy_(&consecutive, &d_one, &child_clique[row + nc * col],
                       &i_one, &clique[i + ldc * j], &i_one);
                break;

              case PackType::Hybrid2: {
                const int nb = S_.blockSize();
                const int jblock_c = col / nb;
                const int jb_c = std::min(nb, nc - nb * jblock_c);
                const int row_ = row - jblock_c * nb;
                const int col_ = col - jblock_c * nb;
                const int start_block_c =
                    clique_block_start_[child_sn][jblock_c];

                const int jblock = j / nb;
                const int jb = std::min(nb, ldc - nb * jblock);
                const int i_ = i - jblock * nb;
                const int j_ = j - jblock * nb;
                const int start_block = clique_block_start_[sn][jblock];

                daxpy_(&consecutive, &d_one,
                       &child_clique[start_block_c + col_ + jb_c * row_], &jb_c,
                       &clique[start_block + j_ + jb * i_], &jb);
              } break;

              case PackType::Hybrid: {
                const int nb = S_.blockSize();
                const int jblock_c = col / nb;
                const int jb_c = std::min(nb, nc - nb * jblock_c);
                const int row_ = row - jblock_c * nb;
                const int col_ = col - jblock_c * nb;
                const int start_block_c =
                    clique_block_start_[child_sn][jblock_c];
                const int ld_c = nc - nb * jblock_c;

                const int jblock = j / nb;
                const int jb = std::min(nb, ldc - nb * jblock);
                const int i_ = i - jblock * nb;
                const int j_ = j - jblock * nb;
                const int start_block = clique_block_start_[sn][jblock];
                const int ld = ldc - nb * jblock;

                daxpy_(&consecutive, &d_one,
                       &child_clique[start_block_c + row_ + ld_c * col_],
                       &i_one, &clique[start_block + i_ + ld * j_], &i_one);
              } break;
            }
            row += consecutive;
          }
        }

        // j < sn_size was already done before, because it was needed before the
        // partial factorisation. Assembling into the clique instead can be done
        // after.
      }
    } else {
      // assemble the child clique into the current clique by blocks of columns.
      // within a block, assemble by rows.

      const int nb = S_.blockSize();
      const int n_blocks = (nc - 1) / nb + 1;

      int row_start{};

      // go through the blocks of columns of the child sn
      for (int b = 0; b < n_blocks; ++b) {
        const int b_start = clique_block_start_[child_sn][b];

        const int col_start = row_start;
        const int col_end = std::min(col_start + nb, nc);

        // go through the rows within this block
        for (int row = row_start; row < nc; ++row) {
          const int i = S_.relindClique(child_sn, row) - sn_size;

          // already assembled into frontal
          if (i < 0) continue;

          // go through the columns of the block
          int col = col_start;
          while (col < col_end) {
            int j = S_.relindClique(child_sn, col);
            if (j < sn_size) {
              ++col;
              continue;
            }
            j -= sn_size;

            // information and sizes of child sn
            const int jblock_c = b;
            const int jb_c = std::min(nb, nc - nb * jblock_c);
            const int row_ = row - jblock_c * nb;
            const int col_ = col - jblock_c * nb;
            const int start_block_c = b_start;

            // sun consecutive entries in a row.
            // consecutive need to be reduced, to account for edge of the block
            const int zeros_stored_row =
                std::max(0, jb_c - (row - row_start) - 1);
            int consecutive = S_.consecutiveSums(child_sn, col);
            const int left_in_child = col_end - col - zeros_stored_row;
            consecutive = std::min(consecutive, left_in_child);

            // consecutive need to account also for edge of block in parent
            const int block_in_parent = j / nb;
            const int col_end_parent =
                std::min((block_in_parent + 1) * nb, ldc);
            const int left_in_parent = col_end_parent - j;
            consecutive = std::min(consecutive, left_in_parent);

            // needed to deal with zeros stored in upper right part of block
            if (consecutive == 0) break;

            // information and sizes of current sn
            const int jblock = block_in_parent;
            const int jb = std::min(nb, ldc - nb * jblock);
            const int i_ = i - jblock * nb;
            const int j_ = j - jblock * nb;
            const int start_block = clique_block_start_[sn][jblock];

            const double d_one = 1.0;
            const int i_one = 1;
            daxpy_(&consecutive, &d_one,
                   &child_clique[start_block_c + col_ + jb_c * row_], &i_one,
                   &clique[start_block + j_ + jb * i_], &i_one);

            col += consecutive;
          }
        }

        row_start += nb;
      }
    }

    // Schur contribution of the child is no longer needed
    delete[] child_clique;

    // move on to the next child
    child_sn = next_children_[child_sn];
  }
  time_assemble_children_C_ += clock.stop();

  return kRetOk;
}

bool Factorise::check() const {
  // Check that the numerical factorisation is correct, by using dense linear
  // algebra operations.
  // Return true if check is successful, or if matrix is too large.
  // To be used for debug.

  if (S_.type() == FactType::AugSys || S_.packFormat() == PackType::Hybrid) {
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

void Factorise::printTimes() const {
  printf("\n----------------------------------------------------\n");
  printf("\t\tFactorise\n");
  printf("----------------------------------------------------\n");
  printf("\nFactorise time          \t%8.4f\n", time_total_);
  printf("\tPrepare:                %8.4f (%4.1f%%)\n", time_prepare_,
         time_prepare_ / time_total_ * 100);
  printf("\tAssembly original:      %8.4f (%4.1f%%)\n", time_assemble_original_,
         time_assemble_original_ / time_total_ * 100);
  printf("\tAssembly into frontal:  %8.4f (%4.1f%%)\n",
         time_assemble_children_F_,
         time_assemble_children_F_ / time_total_ * 100);
  printf("\tAssembly into clique:   %8.4f (%4.1f%%)\n",
         time_assemble_children_C_,
         time_assemble_children_C_ / time_total_ * 100);
  printf("\tDense factorisation:    %8.4f (%4.1f%%)\n", time_factorise_,
         time_factorise_ / time_total_ * 100);

  if (times_dense_fact_[t_dtrsm] + times_dense_fact_[t_dsyrk] +
          times_dense_fact_[t_dgemm] + times_dense_fact_[t_fact] +
          times_dense_fact_[t_dcopy] + times_dense_fact_[t_dscal] +
          times_dense_fact_[t_convert] ==
      0.0) {
    return;
  }

  printf("\t\t  |\n");
  printf("\t\t  |   trsm:     %8.4f (%4.1f%%)\n", times_dense_fact_[t_dtrsm],
         times_dense_fact_[t_dtrsm] / time_factorise_ * 100);
  printf("\t\t  |_  syrk:     %8.4f (%4.1f%%)\n", times_dense_fact_[t_dsyrk],
         times_dense_fact_[t_dsyrk] / time_factorise_ * 100);
  printf("\t\t      gemm:     %8.4f (%4.1f%%)\n", times_dense_fact_[t_dgemm],
         times_dense_fact_[t_dgemm] / time_factorise_ * 100);
  printf("\t\t      fact:     %8.4f (%4.1f%%)\n", times_dense_fact_[t_fact],
         times_dense_fact_[t_fact] / time_factorise_ * 100);
  printf("\t\t      copy:     %8.4f (%4.1f%%)\n", times_dense_fact_[t_dcopy],
         times_dense_fact_[t_dcopy] / time_factorise_ * 100);
  printf("\t\t      copy sch: %8.4f (%4.1f%%)\n",
         times_dense_fact_[t_dcopy_schur],
         times_dense_fact_[t_dcopy_schur] / time_factorise_ * 100);
  printf("\t\t      scal:     %8.4f (%4.1f%%)\n", times_dense_fact_[t_dscal],
         times_dense_fact_[t_dscal] / time_factorise_ * 100);
  printf("\t\t      convert:  %8.4f (%4.1f%%)\n", times_dense_fact_[t_convert],
         times_dense_fact_[t_convert] / time_factorise_ * 100);
}

int Factorise::run(Numeric& num) {
  Clock clock;
  clock.start();

  time_per_sn_.resize(S_.sn());
  times_dense_fact_.resize(TimesInd::t_size);
  clique_block_start_.resize(S_.sn());

  Clock clock_sn;

  int status{};
  for (int sn = 0; sn < S_.sn(); ++sn) {
    clock_sn.start();
    status = processSupernode(sn);
    time_per_sn_[sn] = clock_sn.stop();
    if (status) break;
  }

  time_total_ = clock.stop();

  printTimes();

  if (status) return status;

  // move factorisation to numerical object
  num.sn_columns_ = std::move(sn_columns_);
  num.S_ = &S_;

  return kRetOk;
}