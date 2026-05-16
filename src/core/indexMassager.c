#include <assert.h>

#include "indexMassager.h"
#include "fmpz_mat_helper.h"
#include "fmpz.h"
#include "flint.h"
#include "fmpz_mat.h"
#include "timer.h"
#include "fmpz_rand.h"

/* ceil(log2(x)) for x >= 1; returns 0 for x <= 1. */
static int ceil_log2(int x) {
  int l = 0, v = 1;
  while (v < x) { v <<= 1; l++; }
  return l;
}

/* Compute an index-(m, r) Smith massager (S, U, M, T) for the 2n x 2n
 * matrix B by random projection P = Rem(s * B^{-1} * J, s) followed by
 * a modular Smith form computation.
 *
 * If Q is non-empty, the precomputed projection from
 * largestInvariantFactor is reused (saves one integrality certificate
 * call); otherwise a fresh random J is generated and certified.
 *
 * `kk` is the number of extra projection columns used when reusing Q;
 * otherwise a local k is chosen.
 *
 * Note: the random matrix J is generated with fmpz_mat_randbits
 * (uniform random).  fmpz_mat_randtest must not be used here -- it
 * produces sparse / structured matrices that violate the algorithm's
 * randomness assumption.
 *
 * Returns 1 on success, 0 if a randomized step fails (caller should
 * retry the whole Smith massager attempt). */
int indexMassager(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t M, fmpz_mat_t T, fmpz_mat_t B, int n, int m, int r, fmpz_t s, int kk, fmpz_mat_t Q) {
  assert(m + r <= n);
  if (fmpz_cmp_ui(s, 1) == 0) {
    fmpz_mat_zero(U);
    fmpz_mat_zero(M);
    for (int i = 0; i < T->r; ++i) {
      fmpz_set_ui(fmpz_mat_entry(T, i, i), 1);
      fmpz_set_ui(fmpz_mat_entry(S, i, 0), 1);
    }
    return 1;
  }

  int k = ceil_log2(r) + 4, success = 1;
  fmpz_mat_t J, J1, J2, P, P1, P2;

  fmpz_mat_init(P, 2*n, r+k);
  fmpz_mat_window_init(P1, P, 0, 0, n, r+k);
  fmpz_mat_window_init(P2, P, n, 0, 2*n, r+k);

  if (Q->r == 0) {
    fmpz_mat_init(J, 2*n, r+k);
    fmpz_mat_window_init(J1, J, 0, 0, n, r+k);
    fmpz_mat_window_init(J2, J, n, 0, 2*n, r+k);

    fmpz_mat_zero(J2);
    fmpz_mat_randbits(J1, rand, fmpz_bits(s));
    fmpz_mat_mod(J1, J1, s);
    REAL_TIMER("specialIntCert", success = specialIntCert(P, s, B, J, n, r+k, m));

    if (!fmpz_mat_equal(P2, J2)) {
      printf("Warning: P2 is non-zero\n");
      success = 0;
    }

    fmpz_mat_window_clear(J2);
    fmpz_mat_window_clear(J1);
    fmpz_mat_clear(J);

    if (!success) {
      printf("Warning: sB^{-1}J is not integral Or P2 is non-zero \n");
      goto clear;
    }
  } else {
    k = kk;
    fmpz_mat_window_clear(P1);
    fmpz_mat_window_init(P1, Q, 0, 0, Q->r, Q->c);
  }
  
  REAL_TIMER("computeProjBasis", success = computeProjBasis(S, U, M, T, P1, n, r, s, k));

  fmpz_mat_cmod(M, M, S);
clear:
  fmpz_mat_window_clear(P2);
  fmpz_mat_window_clear(P1);
  fmpz_mat_clear(P);

  return success;
}
