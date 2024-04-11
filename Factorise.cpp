#include "Factorise.h"

#include <fstream>

Factorise::Factorise(const Symbolic& S_input, const int* rowsA_input,
                     const int* ptrA_input, const double* valA_input,
                     int n_input, int nz_input)
    : S{S_input} {
  // Input the symmetric matrix to be factirized in CSC format and the symbolic
  // factorisation coming from Analyze.
  // Only the lower triangular part of the matrix is used.

  if (n_input != S.Size()) {
    printf(
        "Matrix provided to Factorise has size incompatible with symbolic "
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
  ChildrenLinkedList(S.Sn_parent(), firstChildren, nextChildren);

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

void Factorise::ProcessSupernode(int sn) {
  // Assemble frontal matrix for supernode sn, perform partial factorisation and
  // store the result.
  Clock clock;

  clock.start();
  // ===================================================
  // Supernode information
  // ===================================================
  // first and last+1 column of the supernodes
  int sn_begin = S.Sn_start(sn);
  int sn_end = S.Sn_start(sn + 1);
  int sn_size = sn_end - sn_begin;

  // leading dimension of the frontal matrix
  int ldf = S.Ptr(sn + 1) - S.Ptr(sn);

  // leading dimension of the clique matrix
  int ldc = ldf - sn_size;

  // Allocate space for frontal matrix:
  // The front size is ldf and the supernode size is sn_size.
  // The frontal matrix is stored as two dense matrices:
  // Frontal is ldf x sn_size and stores the first sn_size columns that will
  // undergo Cholesky elimination.
  // Clique is ldc x ldc and stores the remaining (ldf - sn_size) columns
  // (without the top part), that do not undergo Cholesky elimination.
  std::vector<double>& Frontal = SnColumns[sn];
  double*& Clique = SchurContribution[sn];

  // Frontal is initialized to zero
  Frontal.resize(ldf * sn_size, 0.0);

  // Clique need not be initialized to zero, provided that the assembly is done
  // properly
  if (ldc > 0) Clique = new double[ldc * ldc];

  time_prepare += clock.stop();

  clock.start();
  // ===================================================
  // Assemble original matrix A into Frontal
  // ===================================================
  // j is relative column index in the frontal matrix
  for (int j = 0; j < sn_size; ++j) {
    // column index in the original matrix
    int col = sn_begin + j;

    // go through the column
    for (int el = ptrA[col]; el < ptrA[col + 1]; ++el) {
      // relative row index in the frontal matrix
      int i = S.Relind_cols(el);

      Frontal[i + j * ldf] = valA[el];
    }
  }
  time_assemble_original += clock.stop();

  // ===================================================
  // Assemble frontal matrices of children into Frontal
  // ===================================================
  clock.start();
  int child_sn = firstChildren[sn];
  while (child_sn != -1) {
    // Schur contribution of the current child
    double*& C = SchurContribution[child_sn];
    if (!C) {
      printf("Error with child supernode\n");
      return;
    }

    // determine size of clique of child
    int child_begin = S.Sn_start(child_sn);
    int child_end = S.Sn_start(child_sn + 1);

    // number of nodes in child sn
    int child_size = child_end - child_begin;

    // size of clique of child sn
    int nc = S.Ptr(child_sn + 1) - S.Ptr(child_sn) - child_size;

    // go through the columns of the contribution of the child
    for (int col = 0; col < nc; ++col) {
      // relative index of column in the frontal matrix
      int j = S.Relind_clique(child_sn, col);

      if (j < sn_size) {
        // assemble into Columns

        // go through the rows of the contribution of the child
        int row = col;
        while (row < nc) {
          // relative index of the entry in the matrix Frontal
          int i = S.Relind_clique(child_sn, row);

          // how many entries to sum
          int consecutive = S.ConsecutiveSums(child_sn, row);

          // use daxpy for summing consecutive entries
          int iOne = 1;
          double dOne = 1.0;
          daxpy(&consecutive, &dOne, &C[row + nc * col], &iOne,
                &Frontal[i + ldf * j], &iOne);
          row += consecutive;
        }
      }

      // If j >= sn_size, we would assemble into Clique.
      // This is delayed until after the partial factorisation, to avoid having
      // to initialize Clique to zero.
    }

    // move on to the next child
    child_sn = nextChildren[child_sn];
  }
  time_assemble_children_F += clock.stop();

  // ===================================================
  // Partial factorisation
  // ===================================================
  clock.start();
  PartialFact_pos_large(ldf, sn_size, Frontal.data(), ldf, Clique, ldc);
  time_factorise += clock.stop();

  // ===================================================
  // Assemble frontal matrices of children into Clique
  // ===================================================
  clock.start();
  child_sn = firstChildren[sn];
  while (child_sn != -1) {
    // Schur contribution of the current child
    double*& C = SchurContribution[child_sn];
    if (!C) {
      printf("Error with child supernode\n");
      return;
    }

    // determine size of clique of child
    int child_begin = S.Sn_start(child_sn);
    int child_end = S.Sn_start(child_sn + 1);

    // number of nodes in child sn
    int child_size = child_end - child_begin;

    // size of clique of child sn
    int nc = S.Ptr(child_sn + 1) - S.Ptr(child_sn) - child_size;

    // go through the columns of the contribution of the child
    for (int col = 0; col < nc; ++col) {
      // relative index of column in the frontal matrix
      int j = S.Relind_clique(child_sn, col);

      if (j >= sn_size) {
        // assemble into Clique

        // adjust relative index to access Clique
        j -= sn_size;

        // go through the rows of the contribution of the child
        int row = col;
        while (row < nc) {
          // relative index of the entry in the matrix Clique
          int i = S.Relind_clique(child_sn, row) - sn_size;

          // how many entries to sum
          int consecutive = S.ConsecutiveSums(child_sn, row);

          // use daxpy for summing consecutive entries
          int iOne = 1;
          double dOne = 1.0;
          daxpy(&consecutive, &dOne, &C[row + nc * col], &iOne,
                &Clique[i + ldc * j], &iOne);
          row += consecutive;
        }
      }

      // j < sn_size was already done before, because it was needed before the
      // partial factorisation. Assembling into the Clique instead can be done
      // after.
    }

    // Schur contribution of the child is no longer needed
    delete[] C;

    // move on to the next child
    child_sn = nextChildren[child_sn];
  }
  time_assemble_children_C += clock.stop();
}

bool Factorise::Check() const {
  // Check that the numerical factorisation is correct, by using dense linear
  // algebra operations.
  // Return true if check is successful, or if matrix is too large.
  // To be used for debug.

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
    for (int col = S.Sn_start(sn); col < S.Sn_start(sn + 1); ++col) {
      // indices to access corresponding entry in the supernode
      int colSn = col - S.Sn_start(sn);
      int ldSn = S.Ptr(sn + 1) - S.Ptr(sn);

      for (int el = S.Ptr(sn); el < S.Ptr(sn + 1); ++el) {
        int row = S.Rows(el);

        // indices to access corresponding entry in the supernode
        int rowSn = el - S.Ptr(sn);

        // skip upper triangle of supernodes
        if (row < col) continue;

        L[row + col * n] = SnColumns[sn][rowSn + colSn * ldSn];
      }
    }
  }

  // Check that sparse factorisation agrees with dense one.
  // This is done by computing the Frobenius norm of the difference between the
  // dense and sparse factors, divided by the Frobenius norm of the dense
  // factor.

  double FrobeniusDense{};
  double FrobeniusDiff{};

  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      double valSparse = L[row + n * col];
      double valDense = M[row + col * n];
      double diff = valSparse - valDense;

      FrobeniusDense += valDense * valDense;
      FrobeniusDiff += diff * diff;
    }
  }

  FrobeniusDense = sqrt(FrobeniusDense);
  FrobeniusDiff = sqrt(FrobeniusDiff);
  check_error = FrobeniusDiff / FrobeniusDense;

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
  printf("\nFactorise time          \t%f\n", time_total);
  printf("\tPrepare:                %f (%4.1f%%)\n", time_prepare,
         time_prepare / time_total * 100);
  printf("\tAssembly original:      %f (%4.1f%%)\n", time_assemble_original,
         time_assemble_original / time_total * 100);
  printf("\tAssembly into Frontal:  %f (%4.1f%%)\n", time_assemble_children_F,
         time_assemble_children_F / time_total * 100);
  printf("\tAssembly into Clique:   %f (%4.1f%%)\n", time_assemble_children_C,
         time_assemble_children_C / time_total * 100);
  printf("\tDense factorisation:    %f (%4.1f%%)\n", time_factorise,
         time_factorise / time_total * 100);
}

void Factorise::Run(Numeric& Num) {
  Clock clock;
  clock.start();

  time_per_Sn.resize(S.Sn());

  Clock clockSn;

  for (int sn = 0; sn < S.Sn(); ++sn) {
    clockSn.start();
    ProcessSupernode(sn);
    time_per_Sn[sn] = clockSn.stop();
  }

  time_total = clock.stop();

  PrintTimes();

  Check();

  // move factorisation to numerical object
  Num.SnColumns = std::move(SnColumns);
  Num.S = &S;
}