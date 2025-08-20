#include <assert.h>

#include "largestInvariantFactor.h"

// Count leading zeros
int clz(unsigned int x)
{
    static const char debruijn32[32] = {
        0, 31, 9, 30, 3, 8, 13, 29, 2, 5, 7, 21, 12, 24, 28, 19,
        1, 10, 4, 14, 6, 22, 25, 20, 11, 15, 23, 26, 16, 27, 17, 18
    };
    x |= x>>1;
    x |= x>>2;
    x |= x>>4;
    x |= x>>8;
    x |= x>>16;
    x++;
    return debruijn32[x*0x076be629>>27];
}

int bits(unsigned int x) {
  return 32 - clz(x);
}

int largestInvariantFactor(fmpz_t s, fmpz_mat_t X, fmpz_mat_t A, fmpz_t startDim) {
  assert(A->c == A->r);

  fmpz_t mx, d;
  fmpz_mat_t B;
  flint_rand_t rand;
  int n = A->c;

  fmpz_init(mx);
  fmpz_init(d);
  fmpz_mat_init(B, n, startDim);
  flint_rand_init(rand);

  fmpz_mat_max(mx, A);
  fmpz_mul_ui(mx, mx, fmpz_bits(mx));
  fmpz_mul_ui(mx, mx, n);
  fmpz_mul_ui(mx, mx, bits(n));
  fmpz_mul_ui(mx, mx, 20);
  fmpz_add_ui(mx, mx, 60);
  fmpz_set_ui(s, 1);
  fmpz_mat_randtest(B, rand, fmpz_bits(mx));

  fmpz_mat_imlSolve(X, d, A, B);
  fmpz_mat_extractLCM(s, X, d);
  fmpz_mat_scalar_mul_fmpz(X, X, s);
  fmpz_mul(s, s, d);
  fmpz_mat_mod(X, X, s);

  flint_rand_clear(rand);
  fmpz_mat_clear(B);
  fmpz_clear(mx);
  fmpz_clear(d);

  return 1;
}
