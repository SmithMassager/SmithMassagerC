#include <assert.h>

#include "indexMassager.h"
#include "fmpz_mat_helper.h"
#include "fmpz.h"
#include "flint.h"
#include "fmpz_mat.h"

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

  int k = 5, success = 1;
  fmpz_mat_t J, J1, J2, P, P1, P2;
  flint_rand_t rand;

  fmpz_mat_init(P, 2*n, r+k);
  fmpz_mat_window_init(P1, P, 0, 0, n, r+k);
  fmpz_mat_window_init(P2, P, n, 0, 2*n, r+k);

  if (Q->r == 0) {
    fmpz_mat_init(J, 2*n, r+k);
    fmpz_mat_window_init(J1, J, 0, 0, n, r+k);
    fmpz_mat_window_init(J2, J, n, 0, 2*n, r+k);
    flint_rand_init(rand);

    fmpz_mat_zero(J2);
    fmpz_mat_randtest(J1, rand, fmpz_bits(s));
    fmpz_mat_mod(J1, J1, s);
    success = specialIntCert(P, s, B, J, n, r+k, m);

    if (!fmpz_mat_equal(P2, J2)) {
      printf("Warning: P2 is non-zero\n");
      success = 0;
    }

    fmpz_mat_window_clear(J2);
    fmpz_mat_window_clear(J1);
    fmpz_mat_clear(J);
    flint_rand_clear(rand);

    if (!success) {
      printf("Warning: sB^{-1}J is not integral Or P2 is non-zero \n");
      goto clear;
    }
  } else {
    k = kk;
    fmpz_mat_window_clear(P1);
    fmpz_mat_window_init(P1, Q, 0, 0, Q->r, Q->c);
  }
  
  success = computeProjBasis(S, U, M, T, P1, n, r, s, k);
  if (!success) {
    printf("computeProjBasis returned false\n");
  }

  fmpz_mat_cmod(M, M, S);
clear:
  fmpz_mat_window_clear(P2);
  fmpz_mat_window_clear(P1);
  fmpz_mat_clear(P);

  return success;
}
