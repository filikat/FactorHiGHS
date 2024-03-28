#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double get_time() {
  struct timespec now;
  clock_gettime(CLOCK_REALTIME, &now);
  return now.tv_sec + now.tv_nsec * 1e-9;
}

#define max(i, j) ((i) >= (j) ? (i) : (j))
#define min(i, j) ((i) >= (j) ? (j) : (i))

// block size
const int nb = 64;

// variables for BLAS calls
double dOne = 1.0;
double dmOne = -1.0;
int iOne = 1;
char LL = 'L';
char NN = 'N';
char RR = 'R';
char TT = 'T';
char UU = 'U';

// BLAS declaration
double ddot(int* n, double* dx, int* incx, double* dy, int* incy);
void dgemv(char* trans, int* m, int* n, double* alpha, double* A, int* lda,
           double* x, int* incx, double* beta, double* y, int* incy);
void dscal(int* n, double* da, double* dx, int* incx);
void dcopy(int* n, double* dx, int* incx, double* dy, int* incy);
void dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha, double* a,
           int* lda, double* beta, double* c, int* ldc);
void dgemm(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
           double* A, int* lda, double* B, int* ldb, double* beta, double* C,
           int* ldc);
void dtrsm(char* side, char* uplo, char* trans, char* diag, int* m, int* n,
           double* alpha, double* a, int* lda, double* b, int* ldb);

// ===========================================================================
// Functions to compute dense partial Cholesky or LDL factorizations with
// left-looking approach, with or without blocking.
//
// A is used to access the first k columns, i.e., M(0:n-1,0:k-1).
// B is used to access the remaining lower triangle, i.e., M(k:n-1,k:n-1).
//
// For indefinite matrices:
// - 2x2 pivoting is not performed. If a zero pivot is found, the code stops.
// - the elements of D are stored as diagonal entries of A; the unit diagonal
//   entries of A are not stored.
//
// These functions are similar to Lapack dpotrf("L", n, A, lda, info)  and
// dpotf2("L", n, A, lda, info), if k >= n.
//
// Arguments:
// - n      : Dimension of matrix M.
// - k      : Number of columns to factorize.
//            If k < n, a partial factorization is computed.
//            If k >= n, a full factorization is computed and B is not used.
// - A      : Array of size (lda * k). To be accessed by columns.
//            On input, it contains the first k columns/rows of M.
//            On output, it contains the trapezoidal factor of the first k
//            columns of M.
// - lda    : Leading dimension of A. It must be at least n. It can be larger
//            if A is stored as a block of a larger matrix.
// - B      : Array of size (ldb * (n-k)). To be accessed by columns.
//            On input, it contains the remaining (n-k) columns of M.
//            On output, it contains the Schur complement.
//            It can be null if k >= n.
// - ldb    : Leading dimension of B. It must be at least (n-k), if k < n. It
//            can be larger if B is stored as a block of a larger martix.
//            Not used if k >= n.
//
// Return:
// - 0      : factorization successful
// - i > 0  : i-th pivot illegal
// - -1     : invalid input
//
// ===========================================================================

int PartialFact_pos_small(int n, int k, double* restrict A, int lda,
                          double* restrict B, int ldb) {
  // ===========================================================================
  // Positive definite factorization without blocks.
  // BLAS calls: ddot, dgemv, dscal, dsyrk.
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("Invalid input to PartialFact\n");
    return -1;
  }

  // quick return
  if (n == 0) return 0;

  // main operations
  for (int j = 0; j < k; ++j) {
    int N = j;
    int M = n - j - 1;

    // update diagonal element
    double Ajj = A[j + lda * j] - ddot(&N, &A[j], &lda, &A[j], &lda);
    if (Ajj <= 0.0 || isnan(Ajj)) {
      A[j + lda * j] = Ajj;
      return j;
    }

    // compute diagonal element
    Ajj = sqrt(Ajj);
    A[j + lda * j] = Ajj;
    double coeff = 1.0 / Ajj;

    // compute column j
    if (j < n - 1) {
      dgemv(&NN, &M, &N, &dmOne, &A[j + 1], &lda, &A[j], &lda, &dOne,
            &A[j + 1 + j * lda], &iOne);
      dscal(&M, &coeff, &A[j + 1 + j * lda], &iOne);
    }
  }

  // update Schur complement
  if (k < n) {
    int N = n - k;
    dsyrk(&LL, &NN, &N, &k, &dmOne, &A[k], &lda, &dOne, B, &ldb);
  }

  return 0;
}

int PartialFact_pos_large(int n, int k, double* restrict A, int lda,
                          double* restrict B, int ldb) {
  // ===========================================================================
  // Positive definite factorization with blocks.
  // BLAS calls: dsyrk, dgemm, dtrsm.
  // ===========================================================================

  double t0, t1, t2, t3, t4, t5, t6, t7;
  double t_copy = 0.0;
  double t_update = 0.0;
  double t_fact = 0.0;
  double t_cols = 0.0;
  double t_schur_copy = 0.0;
  double t_schur = 0.0;

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("Invalid input to PartialFact\n");
    return -1;
  }

  // quick return
  if (n == 0) return 0;

  // j is the starting col of the block of columns
  for (int j = 0; j < k; j += nb) {
    // jb is the size of the block
    int jb = min(nb, k - j);

    // sizes for blas calls
    int N = jb;
    int K = j;
    int M = n - j - jb;

    t0 = get_time();
    // update diagonal block
    dsyrk(&LL, &NN, &N, &K, &dmOne, &A[j], &lda, &dOne, &A[j + lda * j], &lda);
    t1 = get_time();
    // factorize diagonal block
    int info = PartialFact_pos_small(N, N, &A[j + lda * j], lda, NULL, 0);
    if (info != 0) {
      return info + j - 1;
    }
    t2 = get_time();
    if (j + jb < n) {
      // update block of columns
      dgemm(&NN, &TT, &M, &N, &K, &dmOne, &A[j + jb], &lda, &A[j], &lda, &dOne,
            &A[j + jb + lda * j], &lda);

      // solve block of columns with diagonal block
      dtrsm(&RR, &LL, &TT, &NN, &M, &N, &dOne, &A[j + lda * j], &lda,
            &A[j + jb + lda * j], &lda);
    }
    t3 = get_time();
    t_update += t1 - t0;
    t_fact += t2 - t1;
    t_cols += t3 - t2;
  }

  // update Schur complement if partial factorization is required
  if (k < n) {
    int N = n - k;
    t4 = get_time();
    dsyrk(&LL, &NN, &N, &k, &dmOne, &A[k], &lda, &dOne, B, &ldb);
    t5 = get_time();
    t_schur += t5 - t4;
  }

  /*printf("%%%%%%%%%%%%%%%%%%%%%%\n");
  printf("Time profile pos large:\n");
  printf("%15s %12.6f\n", "Time copy", t_copy);
  printf("%15s %12.6f\n", "Time update", t_update);
  printf("%15s %12.6f\n", "Time fact", t_fact);
  printf("%15s %12.6f\n", "Time cols", t_cols);
  printf("%15s %12.6f\n", "Time schur copy", t_schur_copy);
  printf("%15s %12.6f\n", "Time schur", t_schur);
  printf("%%%%%%%%%%%%%%%%%%%%%%\n\n");*/

  return 0;
}

int PartialFact_ind_small(int n, int k, double* restrict A, int lda,
                          double* restrict B, int ldb) {
  // ===========================================================================
  // Infedinite factorization without blocks.
  // BLAS calls: ddot, dgemv, dscal, dcopy, dgemm.
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("Invalid input to PartialFact\n");
    return -1;
  }

  // quick return
  if (n == 0) return 0;

  // main operations
  for (int j = 0; j < k; ++j) {
    int N = j;
    int M = n - j - 1;

    // create temporary copy of row j, multiplied by pivots
    double* temp = malloc(j * sizeof(double));
    for (int i = 0; i < j; ++i) {
      temp[i] = A[j + i * lda] * A[i + i * lda];
    }

    // update diagonal element
    double Ajj = A[j + lda * j] - ddot(&N, &A[j], &lda, temp, &iOne);
    if (Ajj == 0.0 || isnan(Ajj)) {
      A[j + lda * j] = Ajj;
      return j;
    }

    // save diagonal element
    A[j + lda * j] = Ajj;
    double coeff = 1.0 / Ajj;

    // compute column j
    if (j < n - 1) {
      dgemv(&NN, &M, &N, &dmOne, &A[j + 1], &lda, temp, &iOne, &dOne,
            &A[j + 1 + j * lda], &iOne);
      dscal(&M, &coeff, &A[j + 1 + j * lda], &iOne);
    }

    // free temporary copy of row
    free(temp);
  }

  // update Schur complement
  if (k < n) {
    int N = n - k;

    // make temporary copy of M(k:n-1,0:k-1), multiplied by pivots
    double* temp = malloc((n - k) * k * sizeof(double));
    for (int j = 0; j < k; ++j) {
      dcopy(&N, &A[k + j * lda], &iOne, &temp[j * (n - k)], &iOne);
      dscal(&N, &A[j + j * lda], &temp[j * (n - k)], &iOne);
    }

    // update Schur complement using dgemm
    dgemm(&NN, &TT, &N, &N, &k, &dmOne, &A[k], &lda, temp, &N, &dOne, B, &ldb);

    free(temp);
  }

  return 0;
}

int PartialFact_ind_large(int n, int k, double* restrict A, int lda,
                          double* restrict B, int ldb) {
  // ===========================================================================
  // Indefinite factorization with blocks.
  // BLAS calls: dcopy, dscal, dgemm, dtrsm, dsyrk
  // ===========================================================================

  double t0, t1, t2, t3, t4, t5, t6, t7;
  double t_copy = 0.0;
  double t_update = 0.0;
  double t_fact = 0.0;
  double t_cols = 0.0;
  double t_schur_copy = 0.0;
  double t_schur = 0.0;

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("Invalid input to PartialFact\n");
    return -1;
  }

  // quick return
  if (n == 0) return 0;

  // j is the starting col of the block of columns
  for (int j = 0; j < k; j += nb) {
    // jb is the size of the block
    int jb = min(nb, k - j);

    // sizes for blas calls
    int N = jb;
    int K = j;
    int M = n - j - jb;

    t0 = get_time();
    // create temporary copy of block of rows, multiplied by pivots
    double* temp = malloc(j * jb * sizeof(double));
    int ldt = jb;
    for (int i = 0; i < j; ++i) {
      dcopy(&N, &A[j + i * lda], &iOne, &temp[i * ldt], &iOne);
      dscal(&N, &A[i + i * lda], &temp[i * ldt], &iOne);
    }

    t1 = get_time();
    // update diagonal block using dgemm
    dgemm(&NN, &TT, &jb, &jb, &j, &dmOne, &A[j], &lda, temp, &ldt, &dOne,
          &A[j + lda * j], &lda);

    t2 = get_time();
    // factorize diagonal block
    int info = PartialFact_ind_small(N, N, &A[j + lda * j], lda, NULL, 0);
    if (info != 0) {
      return info + j - 1;
    }

    t3 = get_time();
    if (j + jb < n) {
      // update block of columns
      dgemm(&NN, &TT, &M, &N, &K, &dmOne, &A[j + jb], &lda, temp, &ldt, &dOne,
            &A[j + jb + lda * j], &lda);

      // solve block of columns with L
      dtrsm(&RR, &LL, &TT, &UU, &M, &N, &dOne, &A[j + lda * j], &lda,
            &A[j + jb + lda * j], &lda);

      // solve block of columns with D
      for (int i = 0; i < jb; ++i) {
        double coeff = 1.0 / A[j + i + (j + i) * lda];
        dscal(&M, &coeff, &A[j + jb + lda * (j + i)], &iOne);
      }
    }
    t4 = get_time();

    t_copy += t1 - t0;
    t_update += t2 - t1;
    t_fact += t3 - t2;
    t_cols += t4 - t3;

    free(temp);
  }

  // update Schur complement
  if (k < n) {
    int N = n - k;

    t5 = get_time();
    // count number of positive and negative pivots
    int pos_pivot = 0;
    int neg_pivot = 0;
    for (int i = 0; i < k; ++i) {
      if (A[i + lda * i] >= 0.0) {
        ++pos_pivot;
      } else {
        ++neg_pivot;
      }
    }

    // make temporary copies of positive and negative columns separately
    double* temp_pos = malloc((n - k) * pos_pivot * sizeof(double));
    double* temp_neg = malloc((n - k) * neg_pivot * sizeof(double));
    int ldt = n - k;

    // the copies of the columns are multiplied by sqrt(|Ajj|)
    int start_pos = 0;
    int start_neg = 0;
    for (int j = 0; j < k; ++j) {
      double Ajj = A[j + lda * j];
      if (Ajj >= 0.0) {
        Ajj = sqrt(Ajj);
        dcopy(&N, &A[k + j * lda], &iOne, &temp_pos[start_pos * ldt], &iOne);
        dscal(&N, &Ajj, &temp_pos[start_pos * ldt], &iOne);
        ++start_pos;
      } else {
        Ajj = sqrt(-Ajj);
        dcopy(&N, &A[k + j * lda], &iOne, &temp_neg[start_neg * ldt], &iOne);
        dscal(&N, &Ajj, &temp_neg[start_neg * ldt], &iOne);
        ++start_neg;
      }
    }
    t6 = get_time();

    // Update schur complement by subtracting contribution of positive columns
    // and adding contribution of negative columns.
    // In this way, I can use dsyrk instead of dgemm and avoid updating the full
    // square schur complement.
    dsyrk(&LL, &NN, &N, &pos_pivot, &dmOne, temp_pos, &ldt, &dOne, B, &ldb);
    dsyrk(&LL, &NN, &N, &neg_pivot, &dOne, temp_neg, &ldt, &dOne, B, &ldb);
    t7 = get_time();

    t_schur_copy += t6 - t5;
    t_schur += t7 - t6;

    free(temp_pos);
    free(temp_neg);
  }

  printf("%%%%%%%%%%%%%%%%%%%%%%\n");
  printf("Time profile ind large:\n");
  printf("%15s %12.6f\n", "Time copy", t_copy);
  printf("%15s %12.6f\n", "Time update", t_update);
  printf("%15s %12.6f\n", "Time fact", t_fact);
  printf("%15s %12.6f\n", "Time cols", t_cols);
  printf("%15s %12.6f\n", "Time schur copy", t_schur_copy);
  printf("%15s %12.6f\n", "Time schur", t_schur);
  printf("%%%%%%%%%%%%%%%%%%%%%%\n\n");

  return 0;
}
