#pragma once

#include "fmpz_mat.h"
#include "unicert.h"

// Set A = B cmod S.
// A, B have same dimension and S is a column vector of appropiate size.
void fmpz_mat_cmod(fmpz_mat_t A, fmpz_mat_t B, fmpz_mat_t S);
// Set A = B mod p.
void fmpz_mat_mod(fmpz_mat_t A, fmpz_mat_t B, fmpz_t p);
// Set A = B mod p.
void fmpz_mat_mods(fmpz_mat_t A, fmpz_mat_t B, fmpz_t s);
// Set A = |B|.
void fmpz_mat_abs(fmpz_mat_t A, fmpz_mat_t B);
// Set mx = max(|A|).
void fmpz_mat_max(fmpz_t mx, fmpz_mat_t AA);
// imlSolve using fmpz_mat
void fmpz_mat_imlSolve(fmpz_mat_t X, fmpz_t d, fmpz_mat_t A, fmpz_mat_t B);
// Extract the smallest number, a, such that a*X/d is integral.
// Which should be  a = d/gcd(X, d).
void fmpz_mat_extractLCM(fmpz_t a, fmpz_mat_t X, fmpz_t d);
// Set c = max(a, b).
void fmpz_max(fmpz_t c, fmpz_t a, fmpz_t b);
// Calculate res = src^{-1} mod p, where src and res are diagonal matrices represneted by a single column.
void fmpz_mat_modDiagInv(fmpz_mat_t res, fmpz_t p, fmpz_mat_t src);
// Calculate res = A . D mod p where is diagonal.
void fmpz_mat_modDiagMul(fmpz_mat_t res, fmpz_t p, fmpz_mat_t A, fmpz_mat_t D);
// A = (B^T|C^T|D^T)^T
void fmpz_mat_concat_vertical3(fmpz_mat_t A, fmpz_mat_t B, fmpz_mat_t C, fmpz_mat_t D);
int fmpz_uniCert(fmpz_mat_t A);
