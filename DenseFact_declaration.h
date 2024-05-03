#ifndef DENSE_FACT_DECLARATION
#define DENSE_FACT_DECLARATION

#ifdef __cplusplus
extern "C" {
#endif

// dense factorization kernels
int FactPosSmall(char uplo, int n, double* A, int lda);
int FactIndSmall(int n, double* A, int lda);

// dense partial factorization, with blocks
int PartialFactPosLarge(int n, int k, double* A, int lda, double* B, int ldb,
                        double* times);
int PartialFactIndLarge(int n, int k, double* A, int lda, double* B, int ldb,
                        double* times);

// dense partial factorization, in blocked-hybrid format
int PartialFactPosPacked(int n, int k, double* A, int nb, double* B,
                         double* times);

// function to convert A from lower packed, to lower-blocked-hybrid format
void PackedToHybrid(double* A, int nrow, int ncol, int nb);

#ifdef __cplusplus
}
#endif

enum times_ind { t_dtrsm, t_dsyrk, t_dgemm, t_fact, t_dcopy, t_size };

#endif