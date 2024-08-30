#ifndef DENSE_FACT_DECLARATION
#define DENSE_FACT_DECLARATION

#ifdef __cplusplus
extern "C" {
#endif

// dense factorization kernels
int dense_fact_fduf(char uplo, int n, double* A, int lda);
int dense_fact_fiuf(char uplo, int n, double* A, int lda);

// dense partial factorization, with blocks
int dense_fact_pdbf(int n, int k, int nb, double* A, int lda, double* B, int ldb,
                   double* times);
int dense_fact_pibf(int n, int k, int nb, double* A, int lda, double* B, int ldb,
                   double* times);

// dense partial factorization, in blocked-hybrid format
int dense_fact_pdbh(int n, int k, int nb, double* A, double* B, double* times);
int dense_fact_pibh(int n, int k, int nb, double* A, double* B, double* times);

// dense partial factorization, in blocked-hybrid format with hybrid Schur
// complement
int dense_fact_pdbh_2(int n, int k, int nb, double* A, double* B, double* times);
int dense_fact_pibh_2(int n, int k, int nb, double* A, double* B, double* times);

// function to convert A from lower packed, to lower-blocked-hybrid format
int dense_fact_l2h(double* A, int nrow, int ncol, int nb, double* times);

#ifdef __cplusplus
}
#endif

enum TimesInd {
  t_dtrsm,
  t_dsyrk,
  t_dgemm,
  t_fact,
  t_dcopy,
  t_dcopy_schur,
  t_dscal,
  t_convert,
  t_size
};

enum RetValue {
  kRetOk,
  kRetInvalidInput,
  kRetOutOfMemory,
  kRetInvalidPivot,
  kRetGeneric
};

#endif