#ifndef DENSE_FACT_DECLARATION
#define DENSE_FACT_DECLARATION

#ifdef __cplusplus
extern "C" {
#endif

// dense factorization kernels
int dense_fact_fduf(char uplo, int n, double* A, int lda, double thresh,
                    double* regul);
int dense_fact_fiuf(char uplo, int n, double* A, int lda, const int* pivot_sign,
                    double thresh, double* regul, int* n_reg_piv);

// dense partial factorization, with blocks
int dense_fact_pdbf(int n, int k, int nb, double* A, int lda, double* B,
                    int ldb, double thresh, double* regul, double* times);
int dense_fact_pibf(int n, int k, int nb, double* A, int lda, double* B,
                    int ldb, const int* pivot_sign, double thresh,
                    double* regul, int* n_reg_piv, double* times);

// dense partial factorization, in blocked-hybrid format
int dense_fact_pdbh(int n, int k, int nb, double* A, double* B, double thresh,
                    double* regul, double* times);
int dense_fact_pibh(int n, int k, int nb, double* A, double* B,
                    const int* pivot_sign, double thresh, double* regul,
                    int* n_reg_piv, double* times);

// dense partial factorization, in blocked-hybrid format with hybrid Schur
// complement
int dense_fact_pdbs(int n, int k, int nb, double* A, double* B, double thresh,
                    double* regul, double* times);
int dense_fact_pibs(int n, int k, int nb, double* A, double* B,
                    const int* pivot_sign, double thresh, double* regul,
                    int* n_reg_piv, double* times);

// function to convert A from lower packed, to lower-blocked-hybrid format
int dense_fact_l2h(double* A, int nrow, int ncol, int nb, double* times);

#ifdef __cplusplus
}
#endif

#endif