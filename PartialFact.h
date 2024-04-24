#ifndef PARTIALFACT_H
#define PARTIALFACT_H

#ifdef __cplusplus
extern "C" {
#endif

int PartialFactPosLarge(int n, int k, double* A, int lda, double* B, int ldb,
                          double* times);
int PartialFactPosSmall(int n, int k, double* A, int lda, double* B, int ldb);
int PartialFactIndLarge(int n, int k, double* A, int lda, double* B, int ldb);
int PartialFactIndSmall(int n, int k, double* A, int lda, double* B, int ldb);

#ifdef __cplusplus
}
#endif

enum times_ind { t_dtrsm, t_dsyrk, t_dgemm, t_fact, t_size };

#endif