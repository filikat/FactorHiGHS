#ifndef FACTORIZE_H
#define FACTORIZE_H

#include "Auxiliary.h"
#include "PartialFact.h"
#include "Symbolic.h"

class Factorize {
 public:
  // matrix to factorize
  std::vector<int> rowsA{};
  std::vector<int> ptrA{};
  std::vector<double> valA{};
  int n{};
  int nzA{};

  // symbolic factorization
  const Symbolic& S;

  std::vector<double> valL{};

  // children in supernodal elimination tree
  std::vector<int> firstChildren{};
  std::vector<int> nextChildren{};

  // generated elements
  std::vector<double*> SchurContribution{};

  double time_prepare{};
  double time_assemble_original{};
  double time_assemble_children_F{};
  double time_assemble_children_C{};
  double time_factorize{};

 public:
  void Permute(const std::vector<int>& iperm);
  void ProcessSupernode(int sn);
  bool Check() const;

 public:
  Factorize(const Symbolic& S_input, const int* rowsA_input,
            const int* ptrA_input, const double* valA_input, int n_input,
            int nz_input);

  void Run();
};

extern "C" void daxpy(int* n, double* alpha, double* dx, int* incx, double* dy,
                      int* incy);

#endif