#pragma once

#include <inttypes.h>

#include "basis.h"
#include "mpz_matrix.h"

/** \file openblas_residue.h
 * BLAS-compatible matrix
 *
 * An `openblasResidue_t` is a matrix of double-precision floating point numbers
 * representing an integer matrix reduced in the symmetric range with respect
 * to an appropriately-sized modulus. The operations defined here are
 * implemented in terms of OPENBLAS/BLAS routines.  The associated modulus is
 * selected such that the results of underlying BLAS operations can always be
 * represented exactly in the mantissa (53 bits) of a `double` (see basis.h).
 */
struct openblasResidue_t {
  long nrows; /**< row dimension */
  long ncols; /**< column dimension */
  double * _data; /**< data array */
  modulus_t mod; /**< corresponding modulus */
};
typedef struct openblasResidue_t openblasResidue_t;

/** Allocate a new `openblasResidue_t`.
 *
 * The caller is responsible for freeing the new matrix via the `openblas_fini` method.
 * \param nrows row dimension
 * \param ncols column dimension
 * \param p modulus
 */
openblasResidue_t * openblas_init(long nrows, long ncols, long p);

/** Free an `openblasResidue_t` */
void openblas_fini(openblasResidue_t * A);

/* Entry access */
long openblas_getEntry(openblasResidue_t const * A, long idx);
void openblas_setEntry(openblasResidue_t * A, long idx, long val);

/* Setters */
/** Deep matrix copy.
 * Copies the underlying array of doubles from `src` to `dst`.  Source and
 * destination matrices must be of equal dimension.
 */
void openblas_copy(openblasResidue_t * dst, openblasResidue_t const * src);

/** Reduction from arbitrary precision.
 *  Set matrix `dst` to `src` reduced in the symmetric range with respect to `dst`'s modulus.  Matrices are assumed to be of equal dimension.
 */
void openblas_fromMpzMatrix(openblasResidue_t * dst, mpzMatrix_t const * src);

/** Set target to the identity matrix.
 * The target matrix is assumed to be square.
 */
void openblas_identity(openblasResidue_t * A);

/** Zero the target matrix. */
void openblas_zero(openblasResidue_t * A);

/* Query operations*/

/** Check for the zero matrix.
 * \return 1 if matrix A is the zero matrix; 0 otherwise.
 */
int openblas_isZero(openblasResidue_t const * A);

/** Print matrix in a human-readable format.
 */
void openblas_print(FILE * stream, openblasResidue_t const * A);

/* Mutating operations */
/** Scaled matrix addition.
 * Set matrix `dst` to `dst`+ `k`*`src`.
 */
void openblas_add(openblasResidue_t * dst, openblasResidue_t const * src, long k);

/** In-place modular determinant.
 * The input matrix is destroyed in place.
 */
long openblas_determinant(openblasResidue_t * A);

/** In-place modular inverse.
 * Compute the modular inverse of matrix `A` in place.
 * \return
 *  - 1, if the inverse exists
 *  - 0, if the inverse does not exist
 */
int openblas_inverse(openblasResidue_t * A);

/** Matrix multiplication.
 * Set matrix `dst` to `A`*`B`.  Arguments cannot be aliased. Input matrices
 * are assumed to be of compatible dimension.
 */
void openblas_gemm(openblasResidue_t * dst, openblasResidue_t const * A, openblasResidue_t const * B);
void openblas_gemm_helper(openblasResidue_t * dst, openblasResidue_t const * A, openblasResidue_t const * B, long alpha, long beta);

/** Symmetric modular reduction.
 * Reduce (in place) the entries of `dst` in the symmetric range.
 */
void openblas_mods(openblasResidue_t * dst);

/** Convenience quadratic lifting routine.
 * `dst := Xinv.(T - A.M)`
 */
void openblas_quadLift(openblasResidue_t * dst, openblasResidue_t const * T, openblasResidue_t const * A, openblasResidue_t const * M, long xinv);

/** Scalar multiplication.
 * Scale matrix `dst` in place by scalar `k`.
 */
void openblas_scale(openblasResidue_t * dst, long k);

void openblas_set(openblasResidue_t * dst, long r, long c, double src);
double openblas_get(openblasResidue_t * dst, long r, long c);
