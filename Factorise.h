#ifndef FACTORIZE_H
#define FACTORIZE_H

#include "Auxiliary.h"
#include "Numeric.h"
#include "PartialFact.h"
#include "Symbolic.h"

class Factorise {
 public:
  // matrix to factorise
  std::vector<int> rowsA{};
  std::vector<int> ptrA{};
  std::vector<double> valA{};
  int n{};
  int nzA{};

  // symbolic factorisation
  const Symbolic& S;

  // children in supernodal elimination tree
  std::vector<int> firstChildren{};
  std::vector<int> nextChildren{};

  // generated elements
  std::vector<double*> SchurContribution{};

  // columns of L, stored as dense supernodes
  std::vector<std::vector<double>> SnColumns{};

 public:
  void Permute(const std::vector<int>& iperm);
  void ProcessSupernode(int sn);
  bool Check() const;
  void PrintTimes() const;

 public:
  Factorise(const Symbolic& S_input, const int* rowsA_input,
            const int* ptrA_input, const double* valA_input, int n_input,
            int nz_input);

  void Run(Numeric& Num);

  std::vector<double> time_per_Sn{};

  double time_prepare{};
  double time_assemble_original{};
  double time_assemble_children_F{};
  double time_assemble_children_C{};
  double time_factorise{};
  double time_total{};
  mutable double check_error{};
};

extern "C" void daxpy(int* n, double* alpha, double* dx, int* incx, double* dy,
                      int* incy);

#endif