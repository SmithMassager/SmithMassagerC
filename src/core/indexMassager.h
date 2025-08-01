#pragma once

#include "fmpz_mat.h"
#include "fmpz.h"

// IndexMassager(B, n, m, r, s, eps)
//
// Input:
//   B, (2n x 2n) integer matrix with shape:
//         [A              *]
//         [    I_{n-m}     ]
//         [*              *]
//     Last n rows and columns of B^{-1} are integral.
//   n, m, r, positive integer with n+m <= r
//   s, is positive integer multiple of the largest invariant factor of B.
//   eps, in (0,1).
//   Q, an optinal projection matrix.
//
// Output: U, M, T, S
//   U, r x n integer Mod s matrix.
//   M, n x r integer Mod s matrix.
//   T, r x r integer matrix.
//   S, r x r integer matrix, nonsingular and in smith form, with S_rr a divisor
//   of s.
//
//   Index-(m,r) Smith massager (U, M, T, S) for B with T == I_r.
//   With at least 1 - eps probability a maximal index-(m,r) smith massager for
//   B is returned. i.e. S is compromised of the r largest invariant factors of
//   B.
int indexMassager(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t M, fmpz_mat_t T, fmpz_mat_t B, int n, int m, int r, fmpz_t s, int kk, fmpz_mat_t Q);
