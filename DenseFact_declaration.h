#ifndef DENSE_FACT_DECLARATION
#define DENSE_FACT_DECLARATION

#ifdef __cplusplus
extern "C" {
#endif

// dense factorization kernels
int DenseFact_fduf(char uplo, int n, double* A, int lda);
int DenseFact_fiuf(char uplo, int n, double* A, int lda);

// dense partial factorization, with blocks
int DenseFact_pdbf(int n, int k, int nb, double* A, int lda, double* B, int ldb,
                   double* times);
int DenseFact_pibf(int n, int k, int nb, double* A, int lda, double* B, int ldb,
                   double* times);

// dense partial factorization, in blocked-hybrid format
int DenseFact_pdbh(int n, int k, int nb, double* A, double* B, double* times);
int DenseFact_pibh(int n, int k, int nb, double* A, double* B, double* times);

// function to convert A from lower packed, to lower-blocked-hybrid format
int DenseFact_l2h(double* A, int nrow, int ncol, int nb, double* times);

#ifdef __cplusplus
}
#endif

enum times_ind { t_dtrsm, t_dsyrk, t_dgemm, t_fact, t_dcopy, t_dscal, t_size };

enum ret_value {
  ret_ok,
  ret_invalid_input,
  ret_out_of_memory,
  ret_invalid_pivot,
  ret_generic
};

#endif