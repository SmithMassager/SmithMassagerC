#pragma once

#include "fmpz_mat.h"

// SmithMassager(A)
//
// Input:
//  A, n x n integer matrix, nonsingular.
//
// Output: U, M, T, S
//   Tries to find a (0, n)-index Smith Massager for A or false.
//   U, M, T, S, n x n integer matrix.
//   S, in smith form of A.
//   M, with property AM = 0 cmod S
//   T, unit upper triangular with T + UM = 0 cmod S.
int smithMassager(fmpz_mat_t U, fmpz_mat_t M, fmpz_mat_t T, fmpz_mat_t S, fmpz_mat_t A);
int smithMassagerHelper(fmpz_mat_t U, fmpz_mat_t M, fmpz_mat_t T, fmpz_mat_t S, fmpz_mat_t A);
