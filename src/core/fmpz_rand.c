#include "fmpz_rand.h"
#include <time.h>
#include <unistd.h>

int initialized = 0;
flint_rand_t rand;

/* Lazily initialize the global FLINT random state, seeded from
 * time() and getpid() so successive runs follow different
 * randomized paths.  Without this, FLINT uses a fixed default seed
 * and Las Vegas restart counts become deterministic per input. */
void fmpz_rand_init() {
  if (!initialized) {
    flint_rand_init(rand);
    flint_rand_set_seed(rand, (ulong)time(NULL), (ulong)getpid());
    initialized = 1;
  }
}

void fmpz_rand_clear() {
  if (initialized) {
    flint_rand_clear(rand);
  }
}
