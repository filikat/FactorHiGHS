#ifndef GMRES_H
#define GMRES_H

#include <vector>

class LowerMatrix {
  // symmetric matrix in lower triangular form

  const std::vector<int>& rows_;
  const std::vector<int>& ptr_;
  const std::vector<double>& val_;
  const int n;

 public:
  LowerMatrix(const std::vector<int>& rows, const std::vector<int>& ptr,
         const std::vector<double>& val);
  void apply(std::vector<double>& x) const;
};

int Gmres(const LowerMatrix& A, const std::vector<double>& b, std::vector<double>& x,
          double tol, int maxit);

#endif