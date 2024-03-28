#include <fstream>
#include <iostream>
#include <string>

#include "Analyze.h"
#include "Factorize.h"
#include "Highs.h"
#include "io/Filereader.h"

int computeAThetaAT(const HighsSparseMatrix& matrix,
                    const std::vector<double>& theta, HighsSparseMatrix& AAT) {
  // Create a row-wise copy of the matrix
  HighsSparseMatrix AT = matrix;
  AT.ensureRowwise();

  int AAT_dim = matrix.num_row_;
  AAT.num_col_ = AAT_dim;
  AAT.num_row_ = AAT_dim;
  AAT.start_.resize(AAT_dim + 1, 0);

  std::vector<std::tuple<int, int, double>> non_zero_values;

  // First pass to calculate the number of non-zero elements in each column
  //
  int AAT_num_nz = 0;
  std::vector<double> AAT_col_value(AAT_dim, 0);
  std::vector<int> AAT_col_index(AAT_dim);
  std::vector<bool> AAT_col_in_index(AAT_dim, false);
  for (int iRow = 0; iRow < AAT_dim; iRow++) {
    // Go along the row of A, and then down the columns corresponding
    // to its nonzeros
    int num_col_el = 0;
    for (int iRowEl = AT.start_[iRow]; iRowEl < AT.start_[iRow + 1]; iRowEl++) {
      int iCol = AT.index_[iRowEl];
      const double theta_value = !theta.empty() ? theta[iCol] : 1;
      if (!theta_value) continue;
      const double row_value = theta_value * AT.value_[iRowEl];
      for (int iColEl = matrix.start_[iCol]; iColEl < matrix.start_[iCol + 1];
           iColEl++) {
        int iRow1 = matrix.index_[iColEl];
        if (iRow1 < iRow) continue;
        double term = row_value * matrix.value_[iColEl];
        if (!AAT_col_in_index[iRow1]) {
          // This entry is not yet in the list of possible nonzeros
          AAT_col_in_index[iRow1] = true;
          AAT_col_index[num_col_el++] = iRow1;
          AAT_col_value[iRow1] = term;
        } else {
          // This entry is in the list of possible nonzeros
          AAT_col_value[iRow1] += term;
        }
      }
    }
    for (int iEl = 0; iEl < num_col_el; iEl++) {
      int iCol = AAT_col_index[iEl];
      assert(iCol >= iRow);
      const double value = AAT_col_value[iCol];
      if (std::abs(value) > 1e-10) {
        non_zero_values.emplace_back(iRow, iCol, value);
        const int num_new_nz = iRow != iCol ? 2 : 1;
        AAT.start_[iRow + 1]++;
        if (iRow != iCol) AAT.start_[iCol + 1]++;
        AAT_num_nz += num_new_nz;
      }
      AAT_col_value[iCol] =
          0;  // Not strictly necessary, but simplifies debugging
      AAT_col_in_index[iCol] = false;
    }
  }

  // Prefix sum to get the correct column pointers
  for (int i = 0; i < AAT_dim; ++i) AAT.start_[i + 1] += AAT.start_[i];

  AAT.index_.resize(AAT.start_.back());
  AAT.value_.resize(AAT.start_.back());
  AAT.p_end_ = AAT.start_;
  AAT.p_end_.back() = AAT.index_.size();

  std::vector<int> current_positions = AAT.start_;

  // Second pass to actually fill in the indices and values
  for (const auto& val : non_zero_values) {
    int i = std::get<0>(val);
    int j = std::get<1>(val);
    double dot = std::get<2>(val);

    // add dual regularization
    if (i == j) dot += 1;

    AAT.index_[current_positions[i]] = j;
    AAT.value_[current_positions[i]] = dot;
    current_positions[i]++;
    AAT.p_end_[i] = current_positions[i];

    if (i != j) {
      AAT.index_[current_positions[j]] = i;
      AAT.value_[current_positions[j]] = dot;
      current_positions[j]++;
      AAT.p_end_[j] = current_positions[j];
    }
  }
  AAT.p_end_.clear();
  return 0;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Wrong input\n";
    return 1;
  }

  // ===========================================================================
  // Read the problem
  // ===========================================================================

  // Read LP using Highs MPS read
  Highs highs;
  std::string model_file = argv[1];
  HighsStatus status = highs.readModel(model_file);
  assert(status == HighsStatus::kOk);
  status = highs.presolve();
  assert(status == HighsStatus::kOk);
  const HighsLp lp = highs.getPresolvedLp();

  int type = atoi(argv[2]);

  type = 1;

  int nA = lp.a_matrix_.num_col_;
  int mA = lp.a_matrix_.num_row_;
  int nzA = lp.a_matrix_.numNz();

  // ===========================================================================
  // Build the matrix
  // ===========================================================================

  std::vector<int> ptrUpper;
  std::vector<int> rowsUpper;
  std::vector<double> valUpper;
  std::vector<int> ptrLower;
  std::vector<int> rowsLower;
  std::vector<double> valLower;
  int n;
  int nz;

  if (type == 0) {
    // Augmented system

    HighsSparseMatrix At = lp.a_matrix_;
    At.ensureRowwise();

    ptrUpper.resize(nA + mA + 1);
    rowsUpper.resize(nA + nzA + mA);

    // (1,1) diagonal part
    for (int i = 0; i < nA; ++i) {
      ptrUpper[i + 1] = i + 1;
      rowsUpper[i] = i;
    }

    int next = nA;

    // At part
    for (int i = 0; i < mA; ++i) {
      for (int el = At.start_[i]; el < At.start_[i + 1]; ++el) {
        rowsUpper[next] = At.index_[el];
        ++next;
      }
      rowsUpper[next] = nA + i;
      ++next;
      ptrUpper[nA + 1 + i] = next;
    }

    n = nA + mA;
    nz = nA + nzA + mA;
  } else {
    // Normal equations
    std::vector<double> theta;
    HighsSparseMatrix AAt;
    int status = computeAThetaAT(lp.a_matrix_, theta, AAt);

    n = mA;
    nz = AAt.numNz();
    ptrUpper = std::move(AAt.start_);
    rowsUpper = std::move(AAt.index_);
    valUpper = std::move(AAt.value_);
    AAt.clear();
  }

  /*std::vector<int> rowsUpper{0,  0,  1,  2,  3,  4,  3,  4,  5,  3,
                             4,  6,  4,  7,  0,  2,  7,  8,  0,  5,
                             7,  9,  10, 10, 11, 12, 10, 12, 13, 10,
                             12, 14, 2,  7,  12, 15, 10, 12, 14, 16};
  std::vector<int> ptrUpper{0,  1,  3,  4,  5,  6,  9,  12, 14,
                            18, 22, 23, 25, 26, 29, 32, 36, 40};
  std::vector<double> valUpper{20, 1,  20, 30, 35, 20, -1, 1,  25, 1,
                               1,  30, -2, 40, 1,  2,  1,  25, 1,  1,
                               2,  30, 35, 2,  20, 30, 2,  1,  25, 1,
                               1,  20, -1, -1, 2,  40, 1,  2,  -2, 30};
  int n = 17;
  int nz = 40;
  std::vector<int> rowsLower(nz);
  std::vector<int> ptrLower(n + 1);
  std::vector<double> valLower(nz);*/

  // ===========================================================================
  // Symbolic factorization
  // ===========================================================================
  Symbolic S;
  Analyze An(rowsUpper.data(), ptrUpper.data(), n, nz);
  An.Run(S);
  S.Print();

  // ===========================================================================
  // Numerical factorization
  // ===========================================================================

  Factorize F(S, rowsUpper.data(), ptrUpper.data(), valUpper.data(), n, nz);

  F.Run();

  // ===========================================================================
  // Write to file
  // ===========================================================================
  std::ofstream out_file;
  print(out_file, ptrUpper, "ptr");
  print(out_file, rowsUpper, "rows");
  print(out_file, valUpper, "vals");
  print(out_file, S.Perm(), "perm");
  print(out_file, S.Parent(), "parent");
  print(out_file, S.Rowcount(), "rowcount");
  print(out_file, S.Colcount(), "colcount");
  print(out_file, S.Ptr(), "ptrL");
  print(out_file, S.Rows(), "rowsL");
  print(out_file, S.Fsn_start(), "fsn_start");
  print(out_file, F.rowsA, "rowsF");
  print(out_file, F.ptrA, "ptrF");
  print(out_file, F.valA, "valF");

  return 0;
}
