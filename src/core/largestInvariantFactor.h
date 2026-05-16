#pragma once

#include "fmpz_mat.h"
#include "fmpz.h"

/* Compute s, a multiple of the largest invariant factor of A, and a
 * projection X = s * A^{-1} * B mod s for a random B with `startDim`
 * columns.  Returns 1 on success. */
int largestInvariantFactor(fmpz_t s, fmpz_mat_t X, fmpz_mat_t A, slong startDim);
