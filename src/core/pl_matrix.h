#pragma once

#include <inttypes.h>
#include <stdbool.h>

#include "residue.h"
#include "fmpz.h"
#include "fmpz_mat.h"

/** \file pl_matrix.h
 * \brief partial lineraization matrix.
 *
 * A plMatrix_t represents an integer matrix with the associated chosen basis.
 *
 */
struct plMatrix_t {
  long origRows;	  /**<\brief original row dimension */
  long origCols;	  /**<\brief original column dimension */
  long *eparams;           /**<\brief information mapping to the original column and row */
  long X;
  fmpz_mat_t matrix;
};

typedef struct plMatrix_t plMatrix_t;
typedef struct plPLMatrix_t rowPLMatrix_t;
typedef struct plPLMatrix_t colPLMatrix_t;

/* Memory management */
/** Allocate a new `plMatrix_t`.
 *
 * The caller is responsible for freeing the new matrix via the `plMatrix_fini` method.
 * \param orignrows row dimension
 * \param origncols column dimension
 */
plMatrix_t * rowPLMatrix_initFromFmpzMatrix(fmpz_mat_t const src, long X);
plMatrix_t * colPLMatrix_initFromFmpzMatrix(fmpz_mat_t const src, long X);

void plMatrix_fini(plMatrix_t * M);

/* Query operations */
void plMatrix_print(FILE *f, plMatrix_t const * A);
long calculateBase(fmpz_mat_t const src);

/* Mutating operations */
void colCompress(fmpz_mat_t dst, fmpz_mat_t const src, long X, long *eparams, fmpz_mat_struct const *F, bool cmod);
void rowCompress(fmpz_mat_t dst, fmpz_mat_t const src, long X, long *eparams, fmpz_mat_struct const *F, bool cmod);
void MultiplyViaPL(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const B);
void cmodMulviaPL(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const B, fmpz_mat_struct const *F);
void MultiplyViaColPL(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const B);
void cmodMulviaColPL(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const B, fmpz_mat_struct const *F);

void cmod(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const C);
void rmod(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const C);
