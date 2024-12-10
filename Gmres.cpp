#include "Gmres.h"

#include "../ProtoIPM/VectorOperations.h"
#include "Numeric.h"

LowerMatrix::LowerMatrix(const std::vector<int>& rows,
                         const std::vector<int>& ptr,
                         const std::vector<double>& val)
    : rows_{rows}, ptr_{ptr}, val_{val}, n{(int)ptr.size() - 1} {}

void LowerMatrix::apply(std::vector<double>& x) const {
  // Matrix-vector product, with matrix in lower triangular form

  std::vector<double> y(n);

  for (int c = 0; c < n; ++c) {
    for (int el = ptr_[c]; el < ptr_[c + 1]; ++el) {
      const int r = rows_[el];
      const double v = val_[el];

      y[r] += v * x[c];
      if (r != c) y[c] += v * x[r];
    }
  }

  x = std::move(y);
}

void applyRotation(double& x, double& y, double c, double s) {
  double t = c * x + s * y;
  y = -s * x + c * y;
  x = t;
}
void getRotation(double x, double y, double& c, double& s) {
  if (y == 0.0) {
    c = 1.0;
    s = 0.0;
  } else if (std::abs(y) > std::abs(x)) {
    double t = x / y;
    s = 1.0 / sqrt(1.0 + t * t);
    c = t * s;
  } else {
    double t = y / x;
    c = 1.0 / sqrt(1.0 + t * t);
    s = t * c;
  }
}
void update(std::vector<double>& x, int k, std::vector<double>& H, int ldh,
            std::vector<double>& s, std::vector<std::vector<double>>& V) {
  // Solve H * y = s
  std::vector<double> y = s;
  for (int i = k; i >= 0; --i) {
    y[i] /= H[i + ldh * i];
    for (int j = i - 1; j >= 0; --j) {
      y[j] -= H[j + ldh * i] * y[i];
    }
  }
  // x += V * y
  for (int j = 0; j <= k; ++j) vectorAdd(x, V[j], y[j]);
}

int Gmres(const LowerMatrix& A, const Numeric& N, const std::vector<double>& b,
          std::vector<double>& x, double tol, int maxit) {
  // Attempt to solve Ax=b using GMRES, with the preconditioner P given by
  // Numeric object N.
  // Return the number of iterations taken
  //
  // Based on Netlib implementation and Kelley "Iterative Methods for Linear and
  // Nonlinear Equations":
  // https://www.netlib.org/templates/cpp/gmres.h
  // https://www.netlib.org/templates/templates.pdf
  // https://en.wikipedia.org/wiki/Generalized_minimal_residual_method#Example_code

  // sizes
  const int n = x.size();
  const int m = maxit;

  // vectors
  std::vector<std::vector<double>> V;

  // Hessenberg matrix
  std::vector<double> H((m + 1) * m);
  const int ldh = m + 1;

  // Givens rotations
  std::vector<double> sn(m + 1);
  std::vector<double> cs(m + 1);

  // workspace
  std::vector<double> w(n);

  // preconditioned residual P^-1 * (b - A * x)
  std::vector<double> r = x;
  A.apply(r);
  vectorAdd(r, b, 1.0, -1.0);
  N.solve(r);
  double beta = norm2(r);

  // preconditioned rhs P^-1 * b
  w = b;
  N.solve(w);
  double normb = norm2(w);

  if (normb == 0.0) normb = 1.0;
  if (beta / normb <= tol) return 0;

  V.push_back(r);
  vectorScale(V[0], 1.0 / beta);

  // s = beta * e1
  std::vector<double> s(m + 1);
  s[0] = beta;

  // main loop
  int i;
  for (i = 0; i < maxit; ++i) {
    // Arnoldi
    w = V[i];
    A.apply(w);
    N.solve(w);
    V.push_back(w);
    for (int k = 0; k <= i; ++k) {
      H[k + ldh * i] = dotProd(V.back(), V[k]);
      vectorAdd(V.back(), V[k], -H[k + ldh * i]);
    }
    H[i + 1 + ldh * i] = norm2(V.back());

    // re-orthogonalize
    if (norm2(w) + 1e-3 * norm2(V.back()) == norm2(w)) {
      for (int k = 0; k <= i; ++k) {
        double temp = dotProd(V.back(), V[k]);
        H[k + ldh * i] += temp;
        vectorAdd(V.back(), V[k], -temp);
      }
      H[i + 1 + ldh * i] = norm2(V.back());
    }

    // normalize
    vectorScale(V.back(), 1.0 / H[i + 1 + ldh * i]);

    // apply Givens rotations
    for (int k = 0; k < i; ++k)
      applyRotation(H[k + ldh * i], H[k + 1 + ldh * i], cs[k], sn[k]);

    // get latest rotation and apply it also to s
    getRotation(H[i + ldh * i], H[i + 1 + ldh * i], cs[i], sn[i]);
    applyRotation(H[i + ldh * i], H[i + 1 + ldh * i], cs[i], sn[i]);
    applyRotation(s[i], s[i + 1], cs[i], sn[i]);

    printf("%d: %e\n", i, std::abs(s[i + 1]));

    // check termination
    if (std::abs(s[i + 1]) / normb < tol) break;
  }

  // compute solution
  update(x, i, H, ldh, s, V);

  return i + 1;
}
