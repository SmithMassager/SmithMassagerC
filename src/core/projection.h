#pragma once

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_vec.h"
#include "fmpz_mod_mat.h"

// computeProjBasis(P, n, r, s, k)
//
// Input:
//   P, 2n x r+k matrix.
//   n, r, s, positive integer
//
// Pre: The last n rows of P are all 0.
//      k, denotes th extra columns to boost success probability.
//
// Post: S, U, M, T
//   S, r x r integer matrix.
//   U, r x n Mod s matrix.
//   M, n x r Mod s matrix.
//   T, r x r unit upper triangular matrix.
//
//   Let P := [P1], Where P1 and P2 are n x r+k matrices.
//            [P2]
//   -UP1V = F mod s and P1V = MF mod s, where F is in reverse smith form of P1
//   over Mod s and S := sF^-1.
//
//   Calculates the projection basis for Proj(sI, P).
//   Proj(sI, P) := {v length 2n vector | vP/s is integral vector}. And the
//   basis is:
//   [I           ]   [I -M ]    [I    ]
//   [ I          ] * [  I  ]  * [  I  ]
//   [  sF^{-1}   ]   [   I ]    [-U I ]
//   [         I_m]   [    I]    [    I]
int computeProjBasis(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t M, fmpz_mat_t T, fmpz_mat_t P, int n, int r, fmpz_t s, int k);
// SNF is a helper function to do the main work of computeProjBasis.
int SNF(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t M, fmpz_mat_t T, fmpz_mat_t P, fmpz_mat_t E, int n, int r, fmpz_t s, int k, int lower, int upper);
