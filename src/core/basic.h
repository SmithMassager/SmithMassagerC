#pragma once

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "fmpz_vec.h"


void extractGCD(fmpz_t c, fmpz *v, fmpz *w, fmpz_t s, unsigned int n);
void extractMatrixGCDHelper(fmpz *v, fmpz *g, fmpz_t s, fmpz_mat_t A);
void rescale(fmpz_t c, fmpz_t u, fmpz_t g, fmpz_t a, fmpz_t s);
void rescaleE(fmpz_t ret, fmpz_t a, fmpz_t s);
void stab(fmpz_t c, fmpz_t aa, fmpz_t bb, fmpz_t NN);
void vectorGCD(fmpz *v, fmpz_t s, fmpz *w, unsigned int n);

void fmpz_mat_cpy(fmpz_mat_t dst, int dr1, int dc1, int dr2, int dc2, fmpz_mat_t src, int sr1, int sc1, int sr2, int sc2);
void fmpz_mat_multiplyHelper(fmpz_mat_t ret, fmpz_mat_t A, int ar1, int ac1, int ar2, int ac2, fmpz_mat_t B, int br1, int bc1, int br2, int bc2);
void fmpz_mod_mat_multiplyHelper(fmpz_mat_t ret, fmpz_mat_t A, int ar1, int ac1, int ar2, int ac2, fmpz_mat_t B, int br1, int bc1, int br2, int bc2, const fmpz_mod_ctx_t ctx);
void fmpz_mat_scalar_pmod(fmpz_mat_t A, fmpz_mat_t B, fmpz_t mod);
void fmpz_mat_setCol(fmpz_mat_t dst, int col, fmpz *vec, int n);
void fmpz_mat_setRow(fmpz_mat_t dst, int row, fmpz *vec, int n);
void fmpz_mat_vec_vec_mul(fmpz_mat_t dst, int m, int n, fmpz *u, fmpz *v);
void fmpz_max(fmpz_t c, fmpz_t a, fmpz_t b);
