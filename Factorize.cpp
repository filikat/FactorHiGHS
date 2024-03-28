#include "Factorize.h"

#include <fstream>

Factorize::Factorize(const Symbolic& S_input, const int* rowsA_input,
                     const int* ptrA_input, const double* valA_input,
                     int n_input, int nz_input)
    : S{S_input} {
  // Input the symmetric matrix to be factirized in CSC format and the symbolic
  // factorization coming from Analyze.
  // Only the lower triangular part of the matrix is used.
  // Columns are assumed to be sorted.

  if (n_input != S.Size()) {
    printf(
        "Matrix provided to Factorize has size incompatible with symbolic "
        "object.\n");
    return;
  }

  n = n_input;
  rowsA = std::vector<int>(rowsA_input, rowsA_input + nz_input);
  valA = std::vector<double>(valA_input, valA_input + nz_input);
  ptrA = std::vector<int>(ptrA_input, ptrA_input + n_input + 1);

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
  ChildrenLinkedList(S.Fsn_parent(), firstChildren, nextChildren);

  // allocate space for list of generated elements
  SchurContribution.resize(S.Fsn(), {});
}

void Factorize::Permute(const std::vector<int>& iperm) {
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

void Factorize::ProcessSupernode(int sn) {
  // Assemble frontal matrix for supernode sn, perform partial factorization and
  // store the result.

  // =================================================
  // Supernode information
  // =================================================
  // first and last+1 column of the supernodes
  int sn_begin = S.Fsn_start()[sn];
  int sn_end = S.Fsn_start()[sn + 1];
  int sn_size = sn_end - sn_begin;

  // leading dimension of the frontal matrix
  int ldf = S.Ptr()[sn_begin + 1] - S.Ptr()[sn_begin];

  // leading dimension of the clique matrix
  int ldc = ldf - sn_size;

  // Allocate space for frontal matrix:
  // The front size is ldf and the supernode size is sn_size.
  // The frontal matrix is stored as two dense matrices:
  // Frontal is ldf x sn_size and stores the first sn_size columns that will
  // undergo Cholesky elimination.
  // Clique is ldc x ldc and stores the remaining (ldf - sn_size) columns
  // (without the top part), that do not undergo Cholesky elimination.
  std::vector<double> Frontal(ldf * sn_size, 0.0);
  std::vector<double>& Clique = SchurContribution[sn];
  Clique.resize(ldc * ldc, 0.0);

  // =================================================
  // Assemble original matrix A into Frontal
  // =================================================
  // j is relative column index in the frontal matrix
  for (int j = 0; j < sn_size; ++j) {
    // column index in the original matrix
    int col = sn_begin + j;

    // go through the column
    for (int el = ptrA[col]; el < ptrA[col + 1]; ++el) {
      // relative row index in the frontal matrix
      int i = S.Relind_cols()[el];

      Frontal[i + j * ldf] = valA[el];
    }
  }

  // =================================================
  // Assemble frontal matrices of children
  // =================================================
  int child = firstChildren[sn];
  while (child != -1) {
    std::vector<double>& C = SchurContribution[child];
    if (C.empty()) {
      printf("Error with child supernode\n");
      return;
    }
    int nc = sqrt(C.size());

    // go through the columns of the contribution of the child
    for (int col = 0; col < nc; ++col) {
      // relative index of column in the frontal matrix
      int j = S.Relind_clique()[child][col];

      if (j < sn_size) {
        // assemble into Columns

        // go through the rows of the contribution of the child
        for (int row = col; row < nc; ++row) {
          // relative index of the entry in the matrix Frontal
          int i = S.Relind_clique()[child][row];

          Frontal[i + ldf * j] += C[row + nc * col];
        }
      } else {
        // assemble into Clique

        // adjust relative index to access Clique
        j -= sn_size;

        // go through the rows of the contribution of the child
        for (int row = col; row < nc; ++row) {
          // relative index of the entry in the matrix Clique
          int i = S.Relind_clique()[child][row] - sn_size;

          Clique[i + ldc * j] += C[row + nc * col];
        }
      }
    }

    child = nextChildren[child];
  }

  // =================================================
  // Partial factorization
  // =================================================
  PartialFact_pos_large(ldf, sn_size, Frontal.data(), ldf, Clique.data(), ldc);

  for (int j = 0; j < sn_size; ++j) {
    int col = sn_begin + j;
    int i = j;

    for (int el = S.Ptr()[col]; el < S.Ptr()[col + 1]; ++el) {
      valL[el] = Frontal[i + ldf * j];
      ++i;
    }
  }
}

void Factorize::Run() {
  valL.resize(S.Nz());

  for (int sn = 0; sn < S.Fsn(); ++sn) {
    ProcessSupernode(sn);
  }

  std::ofstream out_file;
  print(out_file, valL, "valL");
}