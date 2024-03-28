#ifndef PARTIALFACT_H
#define PARTIALFACT_H

extern "C" int PartialFact_pos_large(int n, int k, double* A, int lda,
                                     double* B, int ldb);
extern "C" int PartialFact_pos_small(int n, int k, double* A, int lda,
                                     double* B, int ldb);
extern "C" int PartialFact_ind_large(int n, int k, double* A, int lda,
                                     double* B, int ldb);
extern "C" int PartialFact_ind_small(int n, int k, double* A, int lda,
                                     double* B, int ldb);

#endif