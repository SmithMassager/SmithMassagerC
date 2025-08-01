#include "fmpz_conv.h"

void fmpzMat_fromMpzMat(fmpz_mat_t Af, mpzMatrix_t *A) {
  for (long i = 0; i < A->nrows; ++i) {
    for (long j = 0; j < A->ncols; ++j) {
      fmpz_set_mpz(fmpz_mat_entry(Af, i, j), mmget(A, i, j));
    }
  }
}

mpzMatrix_t* fmpzMat_toMpzMat(fmpz_mat_t Af) {
  mpzMatrix_t *A;
  A = mpzMatrix_init(fmpz_mat_nrows(Af), fmpz_mat_ncols(Af));
  for (long i = 0; i < A->nrows; ++i) {
    for (long j = 0; j < A->ncols; ++j) {
      fmpz_get_mpz(mmget(A, i, j), fmpz_mat_entry(Af, i, j));
    }
  }
  return A;
}
