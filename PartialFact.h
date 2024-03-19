#ifndef PARTIALFACT_H
#define PARTIALFACT_H

int PartialFact_pos_large(int n, int k, double* A, int lda, double* B, int ldb);
int PartialFact_pos_small(int n, int k, double* A, int lda, double* B, int ldb);
int PartialFact_ind_large(int n, int k, double* A, int lda, double* B, int ldb);
int PartialFact_ind_small(int n, int k, double* A, int lda, double* B, int ldb);

#endif