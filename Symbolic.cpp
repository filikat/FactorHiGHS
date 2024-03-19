#include "Symbolic.h"

#include <iostream>

void Symbolic::Print() const {
  printf("Symbolic factorization:\n");
  printf(" - size %d\n", n);
  printf(" - nonzero entries %d\n", nz);
  printf(" - operations required %.0f\n", operations);
  printf(" - supernodes found %d\n", fsn);
}

int Symbolic::sn_begin(int sn) const {
  // Return the first node in the supernode sn.
  return fsn_start[sn];
}

int Symbolic::sn_end(int sn) const {
  // Return the last node in the supernode sn.
  return fsn_start[sn + 1];
}

int Symbolic::clique_begin(int sn) const {
  // Return the index within rows where the clique corresponding to the
  // supernode sn begins.

  // which column of L to look at
  int j = sn_begin(sn);

  // skip the first entries, because they are nodes in the supernode
  int sn_size = sn_end(sn) - j;
  return ptr[j] + sn_size;
}

int Symbolic::clique_end(int sn) const {
  // Return the index within rows one after the end of the clique of the
  // supernode sn.

  // which column of L to look at
  int j = sn_begin(sn);

  return ptr[j + 1];
}

void Symbolic::clique_info(int sn, int& position, int& snsize,
                           int& cliquesize) const {
  // Given the supernode number sn, position contains the starting entry in rows
  // of the first column of the supernode; snsize is the number of nodes in the
  // supernode; cliquesize is the number of nonzero entries in the remaining of
  // the column.
  // This information corresponds to a node of the clique tree.
  // E.g., after clique_info(4, position, snsize, cliquesize), position is 89,
  // snsize is set to 2 and cliquesize is 3. rows[89:93] contains
  // {8,9,15,17,23}, therefore one can deduce that supernode 4 is made of nodes
  // {8,9} and that the contribution to the Schur complement from supernode 4
  // happens in the rows and columns given by the clique made of nodes
  // {15,17,23}.

  // size of the supernode
  snsize = fsn_start[sn + 1] - fsn_start[sn];

  // first column of the supernode
  int j = fsn_start[sn];

  // number of entries in column j of L, minus snsize
  cliquesize = ptr[j + 1] - ptr[j] - snsize;

  // beginning of column j of L
  position = ptr[j];
}
