#include <stdio.h>
#include <assert.h>

#include "genHerm.h"
#include "fmpz.h"

// return a^b mod p
int power(int a, int b, int p) {
  int ans = 1;
  int x = a;

  while(b) {
    if (b & 1) {
      ans *= x;
      ans %= p;
    }
    x *= x;
    x %= p;
    b = b >> 1;
  }
  return ans;
}

void interestHerm(fmpz_mat_t A) {
  printf("%ld %ld\n", A->c, A->r);
  assert(A->r == A->c);

  fmpz_t p;
  fmpz_init(p);
  fmpz_set_ui(p, A->r);

  assert(fmpz_is_prime(p));
  
  int n = A->r;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      fmpz_set_si(fmpz_mat_entry(A, i, j), power(i, j, n));
    }
  }

  fmpz_clear(p);
}
