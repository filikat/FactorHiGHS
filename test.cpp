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
    if (i == j) dot += 100;

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

  int nA = lp.a_matrix_.num_col_;
  int mA = lp.a_matrix_.num_row_;
  int nzA = lp.a_matrix_.numNz();

  // ===========================================================================
  // Build the matrix
  // ===========================================================================

  std::vector<int> ptrLower;
  std::vector<int> rowsLower;
  std::vector<double> valLower;
  int n;
  int nz;

  if (type == 0) {
    // Augmented system, lower triangular

    const HighsSparseMatrix& A = lp.a_matrix_;

    ptrLower.resize(nA + mA + 1);
    rowsLower.resize(nA + nzA + mA);
    valLower.resize(nA + nzA + mA);

    int next = 0;

    for (int i = 0; i < nA; ++i) {
      // diagonal element
      rowsLower[next] = i;
      valLower[next++] = mA * 1000;

      // column of A
      for (int el = A.start_[i]; el < A.start_[i + 1]; ++el) {
        rowsLower[next] = A.index_[el] + nA;
        valLower[next++] = A.value_[el];
      }

      ptrLower[i + 1] = next;
    }

    // 2,2 block
    for (int i = 0; i < mA; ++i) {
      rowsLower[next] = nA + i;
      valLower[next++] = nA * 1000;
      ptrLower[nA + i + 1] = ptrLower[nA + i] + 1;
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
    ptrLower = std::move(AAt.start_);
    rowsLower = std::move(AAt.index_);
    valLower = std::move(AAt.value_);
    AAt.clear();
  }

  /*int n = 17;
  int nz = 40;
  std::vector<int> rowsLower{0,  1,  8,  9,  1,  2,  8,  15, 3,  5,
                             6,  4,  5,  6,  7,  5,  9,  6,  7,  8,
                             9,  15, 8,  9,  10, 11, 13, 14, 16, 11,
                             12, 13, 14, 15, 16, 13, 14, 16, 15, 16};
  std::vector<int> ptrLower{0,  4,  5,  8,  11, 15, 17, 18, 22,
                            23, 24, 29, 30, 35, 36, 38, 39, 40};
  std::vector<double> valLower{
      20, 1, 1,  1,  20, 20, 1, 1, 20, 1,  1,  20, 1, 1, 1, 20, 1,  20, 20, 1,
      1,  1, 20, 20, 20, 1,  1, 1, 1,  20, 20, 1,  1, 1, 1, 20, 20, 1,  20, 20};*/

  // ===========================================================================
  // Symbolic factorization
  // ===========================================================================
  Symbolic S;
  Analyze An(rowsLower, ptrLower);
  An.Run(S);
  S.Print();

  // ===========================================================================
  // Numerical factorization
  // ===========================================================================

  Factorize F(S, rowsLower.data(), ptrLower.data(), valLower.data(), n, nz);

  F.Run();
  
  // ===========================================================================
  // Write to file
  // ===========================================================================
  std::ofstream out_file;
  print(out_file, ptrLower, "ptr");
  print(out_file, rowsLower, "rows");
  print(out_file, valLower, "vals");
  print(out_file, S.Perm(), "perm");
  print(out_file, S.Sn_start(), "sn_start");
  print(out_file,S.Sn_parent(),"sn_parent");


  return 0;
}
