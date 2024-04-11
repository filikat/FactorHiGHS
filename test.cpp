#include <fstream>
#include <iostream>
#include <random>
#include <regex>
#include <string>

#include "Analyse.h"
#include "Factorise.h"
#include "Highs.h"
#include "hsl_wrapper.h"
#include "io/Filereader.h"

struct MA86Data {
  void* keep;
  ma86_control_d control;
  ma86_info_d info;
  std::vector<int> order;
};

struct MA87Data {
  void* keep;
  ma87_control_d control;
  ma87_info_d info;
  std::vector<int> order;
};

struct MA97Data {
  void* akeep;
  void* fkeep;
  ma97_control_d control;
  ma97_info_d info;
  std::vector<int> order;
};

struct MA57Data {
  void* factors;
  ma57_control_d control;
  ma57_ainfo_d ainfo;
  ma57_finfo_d finfo;
  ma57_sinfo_d sinfo;
  std::vector<int> order;
};

struct MC68Data {
  mc68_control control;
  mc68_info info;
  std::vector<int> order;
};

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
  if (argc < 5) {
    std::cerr << "Wrong input: ./fact pb normEq(0-1) HSL(0-1) Metis(0-1)\n";
    return 1;
  }

  // ===========================================================================
  // Read the problem
  // ===========================================================================

  // Read LP using Highs MPS read
  std::string model_file = argv[1];
  Highs highs;
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
    // Normal equations, full matrix
    std::vector<double> theta;
    HighsSparseMatrix AAt;
    int status = computeAThetaAT(lp.a_matrix_, theta, AAt);

    // extract lower triangle
    n = mA;
    for (int col = 0; col < n; ++col) {
      ptrLower.push_back(valLower.size());
      for (int el = AAt.start_[col]; el < AAt.start_[col + 1]; el++) {
        int row = AAt.index_[el];
        if (row >= col) {
          valLower.push_back(AAt.value_[el]);
          rowsLower.push_back(row);
        }
      }
    }
    ptrLower.push_back(valLower.size());

    nz = ptrLower.back();
  }

  // Test matrix
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
      1,  1, 20, 20, 20, 1,  1, 1, 1,  20, 20, 1,  1, 1, 1, 20, 20, 1,  20,
  20};*/

  // ===========================================================================
  // AMD ordering with MC68
  // ===========================================================================
  std::vector<int> order_to_use{};
  if (atoi(argv[4]) == 0) {
    MC68Data mc68_data;
    mc68_data.order = std::vector<int>(n);
    wrapper_mc68_default_control(&mc68_data.control);
    wrapper_mc68_order(1, n, ptrLower.data(), rowsLower.data(),
                       mc68_data.order.data(), &mc68_data.control,
                       &mc68_data.info);
    if (mc68_data.info.flag < 0) std::cerr << "Error with mc68 ordering\n";
    order_to_use = mc68_data.order;
  }

  // ===========================================================================
  // Symbolic factorisation
  // ===========================================================================
  Symbolic S;
  Analyse An(rowsLower, ptrLower, order_to_use);
  An.Run(S);
  S.Print();

  if (atoi(argv[4]) == 1) order_to_use = An.metis_order;

  // ===========================================================================
  // Numerical factorisation
  // ===========================================================================
  Numeric Num;
  Factorise F(S, rowsLower.data(), ptrLower.data(), valLower.data(), n, nz);
  F.Run(Num);

  // ===========================================================================
  // Solve
  // ===========================================================================

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> distr(-10.0, 10.0);

  // initialize random rhs
  std::vector<double> rhs(n);
  double rhsNorm2{};
  for (int i = 0; i < n; ++i) {
    rhs[i] = distr(rng);
    rhsNorm2 += rhs[i] * rhs[i];
  }
  rhsNorm2 = sqrt(rhsNorm2);

  // solve
  std::vector<double> sol(rhs);
  Num.Solve(sol);

  // ===========================================================================
  // Factorise with MA86
  // ===========================================================================
  double ma86_errorNorm2{};
  double ma86_solNorm2{};
  double ma86_time_analyse{};
  double ma86_time_factorise{};
  if (atoi(argv[3]) == 1) {
    std::vector<double> solMa86(rhs);

    MA86Data ma86_data;

    Clock clock;

    clock.start();
    wrapper_ma86_default_control(&ma86_data.control);
    ma86_data.order = order_to_use;
    wrapper_ma86_analyse(n, ptrLower.data(), rowsLower.data(),
                         ma86_data.order.data(), &ma86_data.keep,
                         &ma86_data.control, &ma86_data.info);
    if (ma86_data.info.flag < 0) std::cerr << "Error with ma86 analyse\n";
    ma86_time_analyse = clock.stop();

    clock.start();
    wrapper_ma86_factor(n, ptrLower.data(), rowsLower.data(), valLower.data(),
                        ma86_data.order.data(), &ma86_data.keep,
                        &ma86_data.control, &ma86_data.info);
    if (ma86_data.info.flag < 0) std::cerr << "Error with ma86 factor\n";
    ma86_time_factorise = clock.stop();

    wrapper_ma86_solve(0, 1, n, solMa86.data(), ma86_data.order.data(),
                       &ma86_data.keep, &ma86_data.control, &ma86_data.info);

    for (int i = 0; i < n; ++i) {
      ma86_errorNorm2 += (solMa86[i] - sol[i]) * (solMa86[i] - sol[i]);
      ma86_solNorm2 += solMa86[i] * solMa86[i];
    }
    ma86_errorNorm2 = sqrt(ma86_errorNorm2);
    ma86_solNorm2 = sqrt(ma86_solNorm2);

    printf("\nRelative error compared to MA86:    %.2e",
           ma86_errorNorm2 / ma86_solNorm2);
    if (ma86_errorNorm2 / ma86_solNorm2 > 1e-6)
      printf(" <================================");
    printf("\n");
    printf("MA86 time analyse:                  %.2e\n", ma86_time_analyse);
    printf("MA86 time factorise:                %.2e\n", ma86_time_factorise);
  }

  // ===========================================================================
  // Factorise with MA87
  // ===========================================================================
  double ma87_errorNorm2{};
  double ma87_solNorm2{};
  double ma87_time_analyse{};
  double ma87_time_factorise{};
  if (atoi(argv[3]) == 1) {
    std::vector<double> solMa87(rhs);

    MA87Data ma87_data;

    Clock clock;

    clock.start();
    wrapper_ma87_default_control(&ma87_data.control);
    ma87_data.order = order_to_use;
    wrapper_ma87_analyse(n, ptrLower.data(), rowsLower.data(),
                         ma87_data.order.data(), &ma87_data.keep,
                         &ma87_data.control, &ma87_data.info);
    if (ma87_data.info.flag < 0) std::cerr << "Error with ma87 analyse\n";
    ma87_time_analyse = clock.stop();

    clock.start();
    wrapper_ma87_factor(n, ptrLower.data(), rowsLower.data(), valLower.data(),
                        ma87_data.order.data(), &ma87_data.keep,
                        &ma87_data.control, &ma87_data.info);
    if (ma87_data.info.flag < 0) std::cerr << "Error with ma87 factor\n";
    ma87_time_factorise = clock.stop();

    wrapper_ma87_solve(0, 1, n, solMa87.data(), ma87_data.order.data(),
                       &ma87_data.keep, &ma87_data.control, &ma87_data.info);

    for (int i = 0; i < n; ++i) {
      ma87_errorNorm2 += (solMa87[i] - sol[i]) * (solMa87[i] - sol[i]);
      ma87_solNorm2 += solMa87[i] * solMa87[i];
    }
    ma87_errorNorm2 = sqrt(ma87_errorNorm2);
    ma87_solNorm2 = sqrt(ma87_solNorm2);

    printf("\nRelative error compared to MA87:    %.2e",
           ma87_errorNorm2 / ma87_solNorm2);
    if (ma87_errorNorm2 / ma87_solNorm2 > 1e-6)
      printf(" <================================");
    printf("\n");
    printf("MA87 time analyse:                  %.2e\n", ma87_time_analyse);
    printf("MA87 time factorise:                %.2e\n", ma87_time_factorise);
  }

  // ===========================================================================
  // Factorise with MA97
  // ===========================================================================
  double ma97_errorNorm2{};
  double ma97_solNorm2{};
  double ma97_time_analyse{};
  double ma97_time_factorise{};
  if (atoi(argv[3]) == 1) {
    std::vector<double> solMa97(rhs);

    MA97Data ma97_data;

    Clock clock;

    clock.start();
    wrapper_ma97_default_control(&ma97_data.control);

    // switch off reordering
    ma97_data.control.ordering = 0;
    ma97_data.order = order_to_use;

    wrapper_ma97_analyse(n, ptrLower.data(), rowsLower.data(),
                         ma97_data.order.data(), &ma97_data.akeep,
                         &ma97_data.control, &ma97_data.info);
    if (ma97_data.info.flag < 0) std::cerr << "Error with ma97 analyse\n";
    ma97_time_analyse = clock.stop();

    clock.start();
    wrapper_ma97_factor(n, ptrLower.data(), rowsLower.data(), valLower.data(),
                        &ma97_data.akeep, &ma97_data.fkeep, &ma97_data.control,
                        &ma97_data.info);
    if (ma97_data.info.flag < 0) std::cerr << "Error with ma97 factor\n";
    ma97_time_factorise = clock.stop();

    wrapper_ma97_solve(0, 1, n, solMa97.data(), &ma97_data.akeep,
                       &ma97_data.fkeep, &ma97_data.control, &ma97_data.info);

    for (int i = 0; i < n; ++i) {
      ma97_errorNorm2 += (solMa97[i] - sol[i]) * (solMa97[i] - sol[i]);
      ma97_solNorm2 += solMa97[i] * solMa97[i];
    }
    ma97_errorNorm2 = sqrt(ma97_errorNorm2);
    ma97_solNorm2 = sqrt(ma97_solNorm2);

    printf("\nRelative error compared to MA97:    %.2e",
           ma97_errorNorm2 / ma97_solNorm2);
    if (ma97_errorNorm2 / ma97_solNorm2 > 1e-6)
      printf(" <================================");
    printf("\n");
    printf("MA97 time analyse:                  %.2e\n", ma97_time_analyse);
    printf("MA97 time factorise:                %.2e\n", ma97_time_factorise);
  }

  // ===========================================================================
  // Factorise with MA57
  // ===========================================================================
  double ma57_errorNorm2{};
  double ma57_solNorm2{};
  double ma57_time_analyse{};
  double ma57_time_factorise{};
  if (atoi(argv[3]) == 1) {
    std::vector<double> solMa57(rhs);

    // ma57 input is in triplet form
    std::vector<int> colsLower(nz);
    for (int col = 0; col < n; ++col) {
      for (int el = ptrLower[col]; el < ptrLower[col + 1]; ++el) {
        colsLower[el] = col;
      }
    }

    MA57Data ma57_data;

    Clock clock;

    clock.start();
    wrapper_ma57_default_control(&ma57_data.control);

    // switch off pivoting
    ma57_data.control.pivoting = 3;

    ma57_data.control.ordering = 2;

    wrapper_ma57_init_factors(&ma57_data.factors);
    wrapper_ma57_analyse(n, nz, rowsLower.data(), colsLower.data(),
                         &ma57_data.factors, &ma57_data.control,
                         &ma57_data.ainfo, ma57_data.order.data());
    if (ma57_data.ainfo.flag < 0) std::cerr << "Error with ma57 analyse\n";
    ma57_time_analyse = clock.stop();

    clock.start();
    wrapper_ma57_factorize(n, nz, rowsLower.data(), colsLower.data(),
                           valLower.data(), &ma57_data.factors,
                           &ma57_data.control, &ma57_data.finfo);
    if (ma57_data.finfo.flag < 0) std::cerr << "Error with ma57 factor\n";
    ma57_time_factorise = clock.stop();

    wrapper_ma57_solve(n, nz, rowsLower.data(), colsLower.data(),
                       valLower.data(), &ma57_data.factors, solMa57.data(),
                       &ma57_data.control, &ma57_data.sinfo);

    for (int i = 0; i < n; ++i) {
      ma57_errorNorm2 += (solMa57[i] - sol[i]) * (solMa57[i] - sol[i]);
      ma57_solNorm2 += solMa57[i] * solMa57[i];
    }
    ma57_errorNorm2 = sqrt(ma57_errorNorm2);
    ma57_solNorm2 = sqrt(ma57_solNorm2);

    printf("\nRelative error compared to MA57:    %.2e",
           ma57_errorNorm2 / ma57_solNorm2);
    if (ma57_errorNorm2 / ma57_solNorm2 > 1e-6)
      printf(" <================================");
    printf("\n");
    printf("MA57 time analyse:                  %.2e\n", ma57_time_analyse);
    printf("MA57 time factorise:                %.2e\n", ma57_time_factorise);
  }

  // ===========================================================================
  // Write to file
  // ===========================================================================

  std::ofstream out_file;
  // print(out_file, ptrLower, "ptr");
  // print(out_file, rowsLower, "rows");
  // print(out_file, valLower, "vals");
  // print(out_file, S.Perm(), "perm");
  print(out_file, S.Sn_start(), "sn_start");
  print(out_file, S.Sn_parent(), "sn_parent");
  print(out_file, S.Ptr(), "ptrsn");
  print(out_file, F.time_per_Sn, "time_per_sn");

  // extract problem name witout mps from path
  std::string pb_name{};
  std::regex rgx("([^/]+)\\.mps");
  std::smatch match;
  std::regex_search(model_file, match, rgx);
  pb_name = match[1];

  // print results
  FILE* file = fopen("results.txt", "a");

  fprintf(
      file, "%15s  |  %12.1e %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f  |  ",
      pb_name.c_str(), An.time_total, An.time_metis / An.time_total * 100,
      An.time_tree / An.time_total * 100, An.time_count / An.time_total * 100,
      An.time_sn / An.time_total * 100, An.time_pattern / An.time_total * 100,
      An.time_relind / An.time_total * 100);

  fprintf(file, "%12.1e %10.1f %10.1f %10.1f %10.1f %10.1f  |  ", F.time_total,
          F.time_prepare / F.time_total * 100,
          F.time_assemble_original / F.time_total * 100,
          F.time_assemble_children_C / F.time_total * 100,
          F.time_assemble_children_F / F.time_total * 100,
          F.time_factorise / F.time_total * 100);

  fprintf(file, "\n");
  fclose(file);

  if (atoi(argv[3])) {
    file = fopen("results_hsl.txt", "a");
    fprintf(file, " %12.5e %12.5e %12.5e %12.5e %12.5e\n", F.time_total,
            ma86_time_factorise, ma87_time_factorise, ma97_time_factorise,
            ma57_time_factorise);
    fclose(file);
  }

  return 0;
}
