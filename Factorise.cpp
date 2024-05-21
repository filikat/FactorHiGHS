#include "Factorise.h"

#include <fstream>

Factorise::Factorise(const Symbolic& S_input,
                     const std::vector<int>& rowsA_input,
                     const std::vector<int>& ptrA_input,
                     const std::vector<double>& valA_input)
    : S{S_input} {
  // Input the symmetric matrix to be factirized in CSC format and the symbolic
  // factorisation coming from Analyze.
  // Only the lower triangular part of the matrix is used.

  n = ptrA_input.size() - 1;

  if (n != S.Size()) {
    printf(
        "Matrix provided to Factorise has size incompatible with symbolic "
        "object.\n");
    return;
  }

  rowsA = rowsA_input;
  valA = valA_input;
  ptrA = ptrA_input;

  // Permute the matrix.
  // This also removes any entry not in the lower triangle.
  Permute(S.Iperm());

  nzA = ptrA.back();

  // Double transpose to sort columns
  std::vector<int> temp_ptr(n + 1);
  std::vector<int> temp_rows(nzA);
  std::vector<double> temp_val(nzA);
  Transpose(ptrA, rowsA, valA, temp_ptr, temp_rows, temp_val);
  Transpose(temp_ptr, temp_rows, temp_val, ptrA, rowsA, valA);

  // create linked lists of children in supernodal elimination tree
  ChildrenLinkedList(S.SnParent(), firstChildren, nextChildren);

  // allocate space for list of generated elements and columns of L
  SchurContribution.resize(S.Sn(), nullptr);
  SnColumns.resize(S.Sn());
}

void Factorise::Permute(const std::vector<int>& iperm) {
  // Symmetric permutation of the lower triangular matrix A based on inverse
  // permutation iperm.
  // The resulting matrix is lower triangular, regardless of the input matrix.

  std::vector<int> work(n, 0);

  // go through the columns to count the nonzeros
  for (int j = 0; j < n; ++j) {
    // get new index of column
    const int col = iperm[j];

    // go through elements of column
    for (int el = ptrA[j]; el < ptrA[j + 1]; ++el) {
      const int i = rowsA[el];

      // ignore potential entries in upper triangular part
      if (i < j) continue;

      // get new index of row
      const int row = iperm[i];

      // since only lower triangular part is used, col is smaller than row
      int actual_col = std::min(row, col);
      ++work[actual_col];
    }
  }

  std::vector<int> new_ptr(n + 1);

  // get column pointers by summing the count of nonzeros in each column.
  // copy column pointers into work
  Counts2Ptr(new_ptr, work);

  std::vector<int> new_rows(new_ptr.back());
  std::vector<double> new_val(new_ptr.back());

  // go through the columns to assign row indices
  for (int j = 0; j < n; ++j) {
    // get new index of column
    const int col = iperm[j];

    // go through elements of column
    for (int el = ptrA[j]; el < ptrA[j + 1]; ++el) {
      const int i = rowsA[el];

      // ignore potential entries in upper triangular part
      if (i < j) continue;

      // get new index of row
      const int row = iperm[i];

      // since only lower triangular part is used, col is smaller than row
      const int actual_col = std::min(row, col);
      const int actual_row = std::max(row, col);

      int pos = work[actual_col]++;
      new_rows[pos] = actual_row;
      new_val[pos] = valA[el];
    }
  }

  ptrA = std::move(new_ptr);
  rowsA = std::move(new_rows);
  valA = std::move(new_val);
}

int Factorise::ProcessSupernode(int sn) {
  // Assemble frontal matrix for supernode sn, perform partial factorisation and
  // store the result.
  Clock clock;

  clock.start();
  // ===================================================
  // Supernode information
  // ===================================================
  // first and last+1 column of the supernodes
  int sn_begin = S.SnStart(sn);
  int sn_end = S.SnStart(sn + 1);
  int sn_size = sn_end - sn_begin;

  // leading dimension of the frontal matrix
  int ldf = S.Ptr(sn + 1) - S.Ptr(sn);

  // leading dimension of the clique matrix
  int ldc = ldf - sn_size;

  // Allocate space for frontal matrix:
  // The front size is ldf and the supernode size is sn_size.
  // The frontal matrix is stored as two dense matrices:
  // frontal is ldf x sn_size and stores the first sn_size columns that will
  // undergo Cholesky elimination.
  // clique is ldc x ldc and stores the remaining (ldf - sn_size) columns
  // (without the top part), that do not undergo Cholesky elimination.
  std::vector<double>& frontal = SnColumns[sn];
  double*& clique = SchurContribution[sn];

  // frontal is initialized to zero
  switch (S.Packed()) {
    case PackType::Full:
    case PackType::Packed:
      frontal.resize(ldf * sn_size, 0.0);
      break;
    case PackType::Hybrid:
    case PackType::Hybrid2:
      frontal.resize(ldf * sn_size - sn_size * (sn_size - 1) / 2, 0.0);
      break;
  }

  // clique need not be initialized to zero, provided that the assembly is done
  // properly
  switch (S.Packed()) {
    case PackType::Full:
      if (ldc > 0) clique = new double[ldc * ldc];
      break;

    case PackType::Packed:
      if (ldc > 0) clique = new double[ldc * (ldc + 1) / 2];
      break;
    case PackType::Hybrid2:
    case PackType::Hybrid: {
      int nb = S.BlockSize();
      int n_blocks = (ldc - 1) / nb + 1;
      clique_block_start[sn].resize(n_blocks + 1);
      int schur_size{};
      for (int j = 0; j < n_blocks; ++j) {
        clique_block_start[sn][j] = schur_size;
        int jb = std::min(nb, ldc - j * nb);
        schur_size += (ldc - j * nb) * jb;
      }
      clique_block_start[sn].back() = schur_size;
      clique = new double[schur_size];
    } break;
  }

  time_prepare += clock.stop();

  clock.start();
  // ===================================================
  // Assemble original matrix A into frontal
  // ===================================================
  // j is relative column index in the frontal matrix
  for (int j = 0; j < sn_size; ++j) {
    // column index in the original matrix
    int col = sn_begin + j;

    // go through the column
    for (int el = ptrA[col]; el < ptrA[col + 1]; ++el) {
      // relative row index in the frontal matrix
      int i = S.RelindCols(el);

      switch (S.Packed()) {
        case PackType::Full:
        case PackType::Packed:
          frontal[i + j * ldf] = valA[el];
          break;
        case PackType::Hybrid:
        case PackType::Hybrid2:
          frontal[i + j * ldf - j * (j + 1) / 2] = valA[el];
          break;
      }
    }
  }
  time_assemble_original += clock.stop();

  // ===================================================
  // Assemble frontal matrices of children into frontal
  // ===================================================
  clock.start();
  int child_sn = firstChildren[sn];
  while (child_sn != -1) {
    // Schur contribution of the current child
    double* child_clique = SchurContribution[child_sn];
    if (!child_clique) {
      printf("Error with child supernode\n");
      return ret_generic;
    }

    // determine size of clique of child
    int child_begin = S.SnStart(child_sn);
    int child_end = S.SnStart(child_sn + 1);

    // number of nodes in child sn
    int child_size = child_end - child_begin;

    // size of clique of child sn
    int nc = S.Ptr(child_sn + 1) - S.Ptr(child_sn) - child_size;

    // go through the columns of the contribution of the child
    for (int col = 0; col < nc; ++col) {
      // relative index of column in the frontal matrix
      int j = S.RelindClique(child_sn, col);

      if (j < sn_size) {
        // assemble into frontal

        // go through the rows of the contribution of the child
        int row = col;
        while (row < nc) {
          // relative index of the entry in the matrix frontal
          int i = S.RelindClique(child_sn, row);

          // how many entries to sum
          int consecutive = S.ConsecutiveSums(child_sn, row);

          // use daxpy for summing consecutive entries
          int i_one = 1;
          double d_one = 1.0;
          switch (S.Packed()) {
            case PackType::Full:
              daxpy(&consecutive, &d_one, &child_clique[row + nc * col], &i_one,
                    &frontal[i + ldf * j], &i_one);
              break;

            case PackType::Packed:
              daxpy(&consecutive, &d_one,
                    &child_clique[row + nc * col - col * (col + 1) / 2], &i_one,
                    &frontal[i + ldf * j], &i_one);
              break;

            case PackType::Hybrid2: {
              int nb = S.BlockSize();
              int jblock = col / nb;
              int jb = std::min(nb, nc - nb * jblock);
              int row_ = row - jblock * nb;
              int col_ = col - jblock * nb;
              int start_block = clique_block_start[child_sn][jblock];
              daxpy(&consecutive, &d_one,
                    &child_clique[start_block + col_ + jb * row_], &jb,
                    &frontal[i + ldf * j - j * (j + 1) / 2], &i_one);
            } break;

            case PackType::Hybrid: {
              int nb = S.BlockSize();
              int jblock = col / nb;
              int row_ = row - jblock * nb;
              int col_ = col - jblock * nb;
              int start_block = clique_block_start[child_sn][jblock];
              int ld = nc - nb * jblock;
              daxpy(&consecutive, &d_one,
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
    child_sn = nextChildren[child_sn];
  }
  time_assemble_children_F += clock.stop();

  // ===================================================
  // Partial factorisation
  // ===================================================
  clock.start();
  switch (S.Packed()) {
    case PackType::Full:
      if (S.Type() == FactType::NormEq) {
        int status = DenseFact_pdbf(ldf, sn_size, S.BlockSize(), frontal.data(),
                                    ldf, clique, ldc, times_dense_fact.data());
        if (status) return status;

      } else {
        int status = DenseFact_pibf(ldf, sn_size, S.BlockSize(), frontal.data(),
                                    ldf, clique, ldc, times_dense_fact.data());
        if (status) return status;
      }
      break;

    case PackType::Packed: {
      // uninitialized temporary full format clique
      double* temp_clique = new double[ldc * ldc];

      if (S.Type() == FactType::NormEq) {
        int status =
            DenseFact_pdbf(ldf, sn_size, S.BlockSize(), frontal.data(), ldf,
                           temp_clique, ldc, times_dense_fact.data());
        if (status) return status;

      } else {
        int status =
            DenseFact_pibf(ldf, sn_size, S.BlockSize(), frontal.data(), ldf,
                           temp_clique, ldc, times_dense_fact.data());
        if (status) return status;
      }

      Clock clock2;
      clock2.start();
      // pack temp_clique into clique
      int pos{};
      for (int j = 0; j < ldc; ++j) {
        int nn = ldc - j;
        int i_one = 1;
        dcopy(&nn, &temp_clique[j + j * ldc], &i_one, &clique[pos], &i_one);
        pos += nn;
      }
      delete[] temp_clique;
      times_dense_fact[t_dcopy] += clock2.stop();
      break;
    }

    case PackType::Hybrid2: {
      int status = DenseFact_l2h(frontal.data(), ldf, sn_size, S.BlockSize(),
                                 times_dense_fact.data());
      if (status) return status;

      if (S.Type() == FactType::NormEq) {
        status = DenseFact_pdbh_2(ldf, sn_size, S.BlockSize(), frontal.data(),
                                  clique, times_dense_fact.data());
        if (status) return status;
      } else {
        status = DenseFact_pibh_2(ldf, sn_size, S.BlockSize(), frontal.data(),
                                  clique, times_dense_fact.data());
        if (status) return status;
      }
    } break;

    case PackType::Hybrid: {
      int status = DenseFact_l2h(frontal.data(), ldf, sn_size, S.BlockSize(),
                                 times_dense_fact.data());
      if (status) return status;

      if (S.Type() == FactType::NormEq) {
        status = DenseFact_pdbh(ldf, sn_size, S.BlockSize(), frontal.data(),
                                clique, times_dense_fact.data());
        if (status) return status;
      } else {
        status = DenseFact_pibh(ldf, sn_size, S.BlockSize(), frontal.data(),
                                clique, times_dense_fact.data());
        if (status) return status;
      }
    } break;
  }

  time_factorise += clock.stop();

  // ===================================================
  // Assemble frontal matrices of children into clique
  // ===================================================
  clock.start();
  child_sn = firstChildren[sn];
  while (child_sn != -1) {
    // Schur contribution of the current child
    double* child_clique = SchurContribution[child_sn];
    if (!child_clique) {
      printf("Error with child supernode\n");
      return ret_generic;
    }

    // determine size of clique of child
    int child_begin = S.SnStart(child_sn);
    int child_end = S.SnStart(child_sn + 1);

    // number of nodes in child sn
    int child_size = child_end - child_begin;

    // size of clique of child sn
    int nc = S.Ptr(child_sn + 1) - S.Ptr(child_sn) - child_size;

    if (S.Packed() != PackType::Hybrid2) {
      //   if (true) {
      //   go through the columns of the contribution of the child
      for (int col = 0; col < nc; ++col) {
        // relative index of column in the frontal matrix
        int j = S.RelindClique(child_sn, col);

        if (j >= sn_size) {
          // assemble into clique

          // adjust relative index to access clique
          j -= sn_size;

          // go through the rows of the contribution of the child
          int row = col;
          while (row < nc) {
            // relative index of the entry in the matrix clique
            int i = S.RelindClique(child_sn, row) - sn_size;

            // how many entries to sum
            int consecutive = S.ConsecutiveSums(child_sn, row);

            // use daxpy for summing consecutive entries
            int i_one = 1;
            double d_one = 1.0;
            switch (S.Packed()) {
              case PackType::Full:
                daxpy(&consecutive, &d_one, &child_clique[row + nc * col],
                      &i_one, &clique[i + ldc * j], &i_one);
                break;

              case PackType::Packed:

              case PackType::Hybrid2: {
                int nb = S.BlockSize();

                int jblock_c = col / nb;
                int jb_c = std::min(nb, nc - nb * jblock_c);
                int row_ = row - jblock_c * nb;
                int col_ = col - jblock_c * nb;
                int start_block_c = clique_block_start[child_sn][jblock_c];

                int jblock = j / nb;
                int jb = std::min(nb, ldc - nb * jblock);
                int i_ = i - jblock * nb;
                int j_ = j - jblock * nb;
                int start_block = clique_block_start[sn][jblock];

                daxpy(&consecutive, &d_one,
                      &child_clique[start_block_c + col_ + jb_c * row_], &jb_c,
                      &clique[start_block + j_ + jb * i_], &jb);
              } break;

              case PackType::Hybrid: {
                int nb = S.BlockSize();

                int jblock_c = col / nb;
                int jb_c = std::min(nb, nc - nb * jblock_c);
                int row_ = row - jblock_c * nb;
                int col_ = col - jblock_c * nb;
                int start_block_c = clique_block_start[child_sn][jblock_c];
                int ld_c = nc - nb * jblock_c;

                int jblock = j / nb;
                int jb = std::min(nb, ldc - nb * jblock);
                int i_ = i - jblock * nb;
                int j_ = j - jblock * nb;
                int start_block = clique_block_start[sn][jblock];
                int ld = ldc - nb * jblock;

                daxpy(&consecutive, &d_one,
                      &child_clique[start_block_c + row_ + ld_c * col_], &i_one,
                      &clique[start_block + i_ + ld * j_], &i_one);
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

      int nb = S.BlockSize();
      int n_blocks = (nc - 1) / nb + 1;

      int row_start{};

      // go through the blocks of columns of the child sn
      for (int b = 0; b < n_blocks; ++b) {
        int b_start = clique_block_start[child_sn][b];

        int col_start = row_start;
        int col_end = std::min(col_start + nb, nc);

        // go through the rows within this block
        for (int row = row_start; row < nc; ++row) {
          int i = S.RelindClique(child_sn, row) - sn_size;

          // already assembled into frontal
          if (i < 0) continue;

          // go through the columns of the block
          int col = col_start;
          while (col < col_end) {
            int j = S.RelindClique(child_sn, col);
            if (j < sn_size) {
              ++col;
              continue;
            }
            j -= sn_size;

            // information and sizes of child sn
            int jblock_c = b;
            int jb_c = std::min(nb, nc - nb * jblock_c);
            int row_ = row - jblock_c * nb;
            int col_ = col - jblock_c * nb;
            int start_block_c = b_start;

            // sun consecutive entries in a row.
            // consecutive need to be reduced, to account for edge of the block
            int zeros_stored_row = std::max(0, jb_c - (row - row_start) - 1);
            int consecutive = S.ConsecutiveSums(child_sn, col);
            int left_in_child = col_end - col - zeros_stored_row;
            consecutive = std::min(consecutive, left_in_child);

            // consecutive need to account also for edge of block in parent
            int block_in_parent = j / nb;
            int col_end_parent = std::min((block_in_parent + 1) * nb, ldc);
            int left_in_parent = col_end_parent - j;
            consecutive = std::min(consecutive, left_in_parent);

            // needed to deal with zeros stored in upper right part of block
            if (consecutive == 0) break;

            // information and sizes of current sn
            int jblock = block_in_parent;
            int jb = std::min(nb, ldc - nb * jblock);
            int i_ = i - jblock * nb;
            int j_ = j - jblock * nb;
            int start_block = clique_block_start[sn][jblock];

            double d_one = 1.0;
            int i_one = 1;
            daxpy(&consecutive, &d_one,
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
    child_sn = nextChildren[child_sn];
  }
  time_assemble_children_C += clock.stop();

  return ret_ok;
}

bool Factorise::Check() const {
  // Check that the numerical factorisation is correct, by using dense linear
  // algebra operations.
  // Return true if check is successful, or if matrix is too large.
  // To be used for debug.

  if (S.Type() == FactType::AugSys || S.Packed() == PackType::Hybrid) {
    printf("\n==> Dense check not available\n");
    return true;
  }

  if (n > 5000) {
    printf("\n==> Matrix is too large for dense check\n\n");
    return true;
  }

  // assemble sparse matrix into dense matrix
  std::vector<double> M(n * n);
  for (int col = 0; col < n; ++col) {
    for (int el = ptrA[col]; el < ptrA[col + 1]; ++el) {
      int row = rowsA[el];

      // insert element in position (row,col)
      M[row + col * n] = valA[el];
    }
  }

  // use Lapack to factorise the dense matrix
  char uplo = 'L';
  int N = n;
  int info;
  dpotrf(&uplo, &N, M.data(), &N, &info);
  if (info != 0) {
    printf("\n==> dpotrf failed\n\n");
    return false;
  }

  // assemble sparse factor into dense factor
  std::vector<double> L(n * n);
  for (int sn = 0; sn < S.Sn(); ++sn) {
    for (int col = S.SnStart(sn); col < S.SnStart(sn + 1); ++col) {
      // indices to access corresponding entry in the supernode
      int col_sn = col - S.SnStart(sn);
      int ldsn = S.Ptr(sn + 1) - S.Ptr(sn);

      for (int el = S.Ptr(sn); el < S.Ptr(sn + 1); ++el) {
        int row = S.Rows(el);

        // indices to access corresponding entry in the supernode
        int row_sn = el - S.Ptr(sn);

        // skip upper triangle of supernodes
        if (row < col) continue;

        L[row + col * n] = SnColumns[sn][row_sn + col_sn * ldsn];
      }
    }
  }

  // Check that sparse factorisation agrees with dense one.
  // This is done by computing the Frobenius norm of the difference between
  // the dense and sparse factors, divided by the Frobenius norm of the dense
  // factor.

  double frobenius_dense{};
  double frobenius_diff{};

  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      double val_sparse = L[row + n * col];
      double val_dense = M[row + col * n];
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

void Factorise::PrintTimes() const {
  printf("\n----------------------------------------------------\n");
  printf("\t\tFactorise\n");
  printf("----------------------------------------------------\n");
  printf("\nFactorise time          \t%8.4f\n", time_total);
  printf("\tPrepare:                %8.4f (%4.1f%%)\n", time_prepare,
         time_prepare / time_total * 100);
  printf("\tAssembly original:      %8.4f (%4.1f%%)\n", time_assemble_original,
         time_assemble_original / time_total * 100);
  printf("\tAssembly into frontal:  %8.4f (%4.1f%%)\n",
         time_assemble_children_F, time_assemble_children_F / time_total * 100);
  printf("\tAssembly into clique:   %8.4f (%4.1f%%)\n",
         time_assemble_children_C, time_assemble_children_C / time_total * 100);
  printf("\tDense factorisation:    %8.4f (%4.1f%%)\n", time_factorise,
         time_factorise / time_total * 100);

  printf("\t\t  |\n");
  printf("\t\t  |   trsm:     %8.4f (%4.1f%%)\n", times_dense_fact[t_dtrsm],
         times_dense_fact[t_dtrsm] / time_factorise * 100);
  printf("\t\t  |_  syrk:     %8.4f (%4.1f%%)\n", times_dense_fact[t_dsyrk],
         times_dense_fact[t_dsyrk] / time_factorise * 100);
  printf("\t\t      gemm:     %8.4f (%4.1f%%)\n", times_dense_fact[t_dgemm],
         times_dense_fact[t_dgemm] / time_factorise * 100);
  printf("\t\t      fact:     %8.4f (%4.1f%%)\n", times_dense_fact[t_fact],
         times_dense_fact[t_fact] / time_factorise * 100);
  printf("\t\t      copy:     %8.4f (%4.1f%%)\n", times_dense_fact[t_dcopy],
         times_dense_fact[t_dcopy] / time_factorise * 100);
  printf("\t\t      copy sch: %8.4f (%4.1f%%)\n",
         times_dense_fact[t_dcopy_schur],
         times_dense_fact[t_dcopy_schur] / time_factorise * 100);
  printf("\t\t      scal:     %8.4f (%4.1f%%)\n", times_dense_fact[t_dscal],
         times_dense_fact[t_dscal] / time_factorise * 100);
  printf("\t\t      convert:  %8.4f (%4.1f%%)\n", times_dense_fact[t_convert],
         times_dense_fact[t_convert] / time_factorise * 100);
}

int Factorise::Run(Numeric& Num) {
  Clock clock;
  clock.start();

  time_per_Sn.resize(S.Sn());
  times_dense_fact.resize(times_ind::t_size);
  clique_block_start.resize(S.Sn());

  Clock clock_sn;

  int status{};
  for (int sn = 0; sn < S.Sn(); ++sn) {
    clock_sn.start();
    status = ProcessSupernode(sn);
    time_per_Sn[sn] = clock_sn.stop();
    if (status) break;
  }

  time_total = clock.stop();

  PrintTimes();

  if (status) return status;

  // move factorisation to numerical object
  Num.SnColumns = std::move(SnColumns);
  Num.S = &S;

  return ret_ok;
}