#pragma once

#include "fmpz_mat.h"
#include "mpz_matrix.h"

// fastIntCertificate(s, A, B, n, m)
//
// Input:
//   A, n x n matrix over Z, nonsingular.
//   s, positive integer.
//   B, n x m matrix over Z(s).
//   n, m, positive integer.
//
// Output:
//   If sA^-1B is integral then return modp(sA^-1B, s) otherwise return false.
int fastIntCert(mpz_t s, mpzMatrix_t *Res, mpzMatrix_t *A, mpzMatrix_t *B, int n, int m);

int fmpzFastIntCert(fmpz_t s, fmpz_mat_t Res, fmpz_mat_t A, fmpz_mat_t B, int n, int m);
