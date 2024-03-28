#include <fstream>
#include <iostream>
#include <string>

#include "Analyze.h"
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
    if (i == j) dot += 0;

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
  int type = 1;

  std::ifstream names("../../Netlib/netlib_names.txt");
  std::string pb;
  double total_time{};

  while (getline(names, pb)) {
    printf("\n====================\n");
    printf("Running %s\n\n", pb.c_str());

    // ===========================================================================
    // Read the problem
    // ===========================================================================

    // Read LP using Highs MPS read
    Highs highs;
    char pb_path[80];
    snprintf(pb_path, 80, "../../Netlib/data/%s", pb.c_str());
    HighsStatus status = highs.readModel(pb_path);
    assert(status == HighsStatus::kOk);
    status = highs.presolve();
    assert(status == HighsStatus::kOk);
    const HighsLp lp = highs.getPresolvedLp();

    int nA = lp.a_matrix_.num_col_;
    int mA = lp.a_matrix_.num_row_;
    int nzA = lp.a_matrix_.numNz();

    // ===========================================================================
    // Build the matrix
    // ===========================================================================

    std::vector<int> ptr;
    std::vector<int> rows;
    int n;
    int nz;

    if (type == 0) {
      // Augmented system

      HighsSparseMatrix At = lp.a_matrix_;
      At.ensureRowwise();

      ptr.resize(nA + mA + 1);
      rows.resize(nA + nzA + mA);

      // (1,1) diagonal part
      for (int i = 0; i < nA; ++i) {
        ptr[i + 1] = i + 1;
        rows[i] = i;
      }

      int next = nA;

      // At part
      for (int i = 0; i < mA; ++i) {
        for (int el = At.start_[i]; el < At.start_[i + 1]; ++el) {
          rows[next] = At.index_[el];
          ++next;
        }
        rows[next] = nA + i;
        ++next;
        ptr[nA + 1 + i] = next;
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
      ptr = std::move(AAt.start_);
      rows = std::move(AAt.index_);
      AAt.clear();
    }

    // ===========================================================================
    // Symbolic factorization
    // ===========================================================================
    Clock clock;
    clock.start();
    Symbolic S;
    Analyze An(rows.data(), ptr.data(), n, nz);
    An.Run(S);
    total_time += clock.stop();
    S.Print();
  }

  printf("\nTotal time %f\n", total_time);

  return 0;
}
