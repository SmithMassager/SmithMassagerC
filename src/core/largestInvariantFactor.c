#include <assert.h>

#include "largestInvariantFactor.h"
#include "timer.h"
#include "basic.h"
#include "fmpz_rand.h"


/* Compute s, a multiple of the largest invariant factor of A, and
 * the projection X = s * A^{-1} * B mod s for a random B with
 * `startDim` columns.
 *
 * Note: the random matrix is generated with fmpz_mat_randbits (uniform
 * random).  fmpz_mat_randtest must not be used here -- it is FLINT's
 * edge-case generator and produces sparse / structured matrices that
 * violate the algorithm's randomness assumption. */
int largestInvariantFactor(fmpz_t s, fmpz_mat_t X, fmpz_mat_t A, slong startDim) {
  assert(A->c == A->r);

  fmpz_t mx, d;
  fmpz_mat_t B;
  int n = A->c;

  fmpz_init(mx);
  fmpz_init(d);
  fmpz_mat_init(B, n, startDim);

  fmpz_mat_max(mx, A);
  fmpz_mul_ui(mx, mx, fmpz_bits(mx));
  fmpz_mul_ui(mx, mx, n);
  fmpz_mul_ui(mx, mx, bits(n));
  fmpz_mul_ui(mx, mx, 20);
  fmpz_add_ui(mx, mx, 60);
  fmpz_set_ui(s, 1);
  fmpz_mat_randbits(B, rand, fmpz_bits(mx));

  REAL_TIMER("ImlSolve in largestInvariantFactor", fmpz_mat_imlSolve(X, d, A, B));
  fmpz_mat_extractLCM(s, X, d);
  fmpz_mat_scalar_mul_fmpz(X, X, s);
  fmpz_mul(s, s, d);
  fmpz_mat_mod(X, X, s);

  fmpz_mat_clear(B);
  fmpz_clear(mx);
  fmpz_clear(d);

  return 1;
}
