#include "fmpz_rand.h"

int initialized = 0;
flint_rand_t rand;

void fmpz_rand_init() {
  if (!initialized) {
    flint_rand_init(rand);
  }
}

void fmpz_rand_clear() {
  if (initialized) {
    flint_rand_clear(rand);
  }
}
