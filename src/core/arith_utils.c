#include "arith_utils.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gmp.h"

  /* Return g = gcd(a,b) and s, t, u, v with
   [ g ] = [ s  t ][ a ]
   [ 0 ]   [ u  v ][ b ],
   v > 0 and 0 <= t < v.
  */
void gcdex(mpz_t g, mpz_t s, mpz_t t, mpz_t u, mpz_t v, mpz_t const a, mpz_t const b)
{
  assert(g != a);
  assert(g != b);

  /* g = sa + tb */
  mpz_gcdext(g, t, s, b, a);

  /* v := |a / g| */
  mpz_tdiv_q(v, a, g);
  mpz_abs(v, v);

  /* u := -vb/a */
  mpz_mul(u, v, b);
  mpz_tdiv_q(u, u, a);
  mpz_neg(u, u);

  /* t := t mod v (enforce 0 <= t < v) */
  mpz_mod(t, t, v);

  /* s := (g - tb) / a (update s to reflect change to t)*/
  mpz_set(s, g);
  mpz_submul(s, t, b);
  mpz_tdiv_q(s, s, a);
}

// minimal c s.t. gcd(a, b, N) = gcd(a + cb, N) 
// Micmic the implementation of Stab as in the maple file.
void stab(mpz_t c, mpz_t const a, mpz_t const b, mpz_t const N) {
  mpz_t g, aa, bb, NN;
  mpz_init(g); mpz_init_set(aa, a); mpz_init_set(bb, b); mpz_init_set(NN, N);
  mpz_set_ui(c, 0);
  
  mpz_gcd(g, aa, bb);
  if (mpz_cmp_ui(g, 0) == 0) { mpz_clears(g, aa, bb, NN, NULL); return; }
  mpz_tdiv_q(aa, aa, g);
  mpz_tdiv_q(bb, bb, g);
  mpz_gcd(g, aa, bb);
  if (mpz_cmp_ui(g, 0) == 0) { mpz_clears(g, aa, bb, NN, NULL); return; }

  mpz_tdiv_q(NN, NN, g);
  mpz_tdiv_q(aa, aa, g);
  mpz_tdiv_q(bb, bb, g);

  unsigned long end = mpz_get_ui(NN);
  for (unsigned long i = 0; i < end; ++i) {
    mpz_gcd(g, aa, NN);
    if (mpz_cmp_ui(g, 1) == 0) { mpz_set_ui(c, i); mpz_clears(g, aa, bb, NN, NULL); return; }
    mpz_add(aa, aa, bb);
  }
}

// Set e s.t. (a * e) = (a, s)
void rescale(mpz_t e, mpz_t const a, mpz_t const s) {
  // Calculate e = stab(u + c * s/g) where g = gcd(a, s), g = au + sv, c = stab(u, s/g, s).
  mpz_t u, v, g, tmp1, tmp2;
  mpz_inits(u, v, g, tmp1, tmp2, NULL);
  mpz_gcdext(g, u, v, a, s);
  mpz_tdiv_q(tmp1, s, g);
  stab(tmp2, u, tmp1, s);
  mpz_mul(tmp2, tmp1, tmp2);
  mpz_add(e, u, tmp2);
  mpz_clears(u, v, g, tmp1, tmp2, NULL);
}

/* Return 1/n mod p */
long modInverseL(long nn, long pp)
{
  long ret;
  mpz_t n;
  mpz_init_set_ui(n, nn);
  ret = modInverseMpz(n, pp);

  mpz_clear(n);
  return ret;
}

long modInverseMpz(mpz_t const n, long pp)
{
  long ret;
  mpz_t rslt, p;
  mpz_init(rslt);
  mpz_init_set_ui(p, pp);
  mpz_invert(rslt, n, p);

  assert(mpz_fits_slong_p(rslt));
  ret = mpz_get_si(rslt);

  mpz_clear(rslt);
  mpz_clear(p);
  return ret;
}

/* Return n mod p in the symmetric range. */
long modsL(long n, long p)
{
  long half_p = p/2;
  long r = n % p;
  if (r < -half_p) {
    r += p;
  } else if (r > half_p) {
    r -= p;
  }
  return r;
}

long modsMpzL(mpz_t const r, long p)
{
  long half_p = p/2;

  long rr = mpz_fdiv_ui(r, p);

  if ( rr > half_p ) {
    rr -= p;
  }
  return rr;
}

void modsMpzMpz(mpz_t r, mpz_t const p)
{
  mpz_t half_p;
  mpz_init(half_p);
  mpz_fdiv_q_2exp(half_p, p, 1);

  mpz_fdiv_r(r, r, p);
  if ( 0 < mpz_cmp(r, half_p)) {
    mpz_sub(r, r, p);
  }
  mpz_clear(half_p);
}

/* Return next (previous) prime. */
long nextprime(long nn)
{
  long ret;
  mpz_t n;
  mpz_init_set_ui(n, nn);
  mpz_nextprime(n, n);
  assert(mpz_fits_ulong_p(n));
  ret = mpz_get_ui(n);

  mpz_clear(n);
  return ret;
}

long prevprime(long nn)
{
  long ret;
  mpz_t n;
  mpz_init_set_ui(n, nn);
  mpz_sub_ui(n, n, 1);
  while(!mpz_probab_prime_p(n, 10)) {
    mpz_sub_ui(n, n, 1);
  }
  assert(mpz_fits_ulong_p(n));
  ret = mpz_get_ui(n);

  mpz_clear(n);
  return ret;
}

void addmul_si(mpz_t ret, mpz_t const a, long b)
{
  if(b >= 0) {
    mpz_addmul_ui(ret, a, b);
  } else {
    mpz_submul_ui(ret, a, -b);
  }
}

// Calculate ret = a + c b mod n
void mpz_addmul_mod(mpz_t ret, mpz_t const a, mpz_t const c, mpz_t const b, mpz_t const n)
{
  mpz_t tmp;
  mpz_init(tmp);
  mpz_mul(tmp, c, b);
  mpz_add(ret, a, tmp);
  mpz_clear(tmp);
}
