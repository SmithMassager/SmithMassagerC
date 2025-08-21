#include <unistd.h>
#include "basic.h"

// fmpz_stab(c, a, b, N)
//
// Input: integers a, b, N, N > 0
//
// Post: nonegative integer c such that
//   gcd(a, b, N) = gcd(a + c*b, N)
//   c is minimal, i.e. c <= ceil(2 * log(N2)^3/2)
void fmpz_stab(fmpz_t c, fmpz_t aa, fmpz_t bb, fmpz_t NN) {
  fmpz_t g, a, b, N;

  fmpz_init(g);
  fmpz_init(a);
  fmpz_init(b);
  fmpz_init(N);

  fmpz_set(a, aa);
  fmpz_set(b, bb);
  fmpz_set(N, NN);

  fmpz_gcd3(g, a, b, N);
  if (fmpz_equal_ui(g, 0)) {
    fmpz_set_ui(c, 0);
    goto clear;
  }

  fmpz_tdiv_q(a, a, g);
  fmpz_tdiv_q(b, b, g);
  fmpz_tdiv_q(N, N, g);
  fmpz_gcd(g, a, b);
  if (fmpz_equal_ui(g, 0)) {
    fmpz_set_ui(c, 0);
    goto clear;
  }

  fmpz_gcd(g, a, N);
  if (fmpz_equal_ui(g, 1)) {
    fmpz_set_ui(c, 0);
    goto clear;
  }
  fmpz_set_ui(c, 1);
  while(1) {
    fmpz_add(a, a, b);
    fmpz_gcd(g, a, N);
    if (fmpz_equal_ui(g,1)) {
      goto clear;
    }
    fmpz_add_ui(c, c, 1);
  }
  
clear:
  // Clear a,b,N;
  fmpz_clear(g);
  fmpz_clear(a);
  fmpz_clear(b);
  fmpz_clear(N);
  return;
}

// vectorGCD(v, s, w, n)
//
// Input:
//   w, nx1 mod s vector
//   s, positive integer
//
// Pre: Assume v and w are appropiate size.
//
// Post: v
//   <v, w> = gcd(w)
void vectorGCD(fmpz *v, fmpz_t s, fmpz *w, unsigned int n) {
  // TODO: I think we should create a fmpz_vec_t type which is just a pointer to *fmpz_vec struct which contains number of entries and a pointer to *fmpz
  fmpz_t g;
  fmpz_t c;
  fmpz_t tmp;

  fmpz_init(g);
  fmpz_init(c);
  fmpz_init(tmp);

  fmpz_mod(g, &w[0], s);
  fmpz_set_ui(&v[0], 1);

  for (int i = 1; i < n; ++i) {
    if (fmpz_equal_ui(&w[i], 0)) continue;
    if (fmpz_equal_ui(g, 1)) break;
    fmpz_stab(c, g, &w[i], s);
    fmpz_set(&v[i], c);
    fmpz_mul(tmp, c, &w[i]);
    fmpz_add(g, g, tmp);
    fmpz_mod(g, g, s);
  }

  fmpz_clear(g);
  fmpz_clear(c);
  fmpz_clear(tmp);
}

// extractGCD(c, v, w, s, n)
//
// Input:
//   v, w, vector Mod s of same size n.
//   s, positive integer.
//
// Post: c
//   s.t. gcd(v + cw, s) = gcd(v, w, s).
void extractGCD(fmpz_t c, fmpz *v, fmpz *w, fmpz_t s, unsigned int n) {
  fmpz_t g;
  fmpz m[2];
  fmpz tmp[2];

  fmpz_init(g);
  fmpz_init(&m[0]);
  fmpz_init(&m[1]);
  fmpz_init(&tmp[0]);
  fmpz_init(&tmp[1]);

  _fmpz_vec_content(g, w, n);
  fmpz_gcd(&m[0], &v[0], s);
  fmpz_gcd(&m[1], &w[0], s);
  for (int i = 1; i < n; ++i) {
    if (fmpz_equal(g, &m[1])) {
      break;
    }
    fmpz_stab(c, &m[1], &w[i], s);
    fmpz_set(&tmp[0], &v[i]);
    fmpz_set(&tmp[1], &w[i]);
    _fmpz_vec_scalar_addmul_fmpz(m, tmp, 2, c);
    _fmpz_vec_scalar_mod_fmpz(m, m, 2, s);
  }

  fmpz_stab(c, &m[0], &m[1], s);

  fmpz_clear(g);
  fmpz_clear(&m[0]);
  fmpz_clear(&m[1]);
  fmpz_clear(&tmp[0]);
  fmpz_clear(&tmp[1]);
}

void fmpz_vec_gcd(fmpz_t g, fmpz *v, int n) {
  fmpz_set(g, &v[0]);
  for (int i = 1; i < n; ++i) {
    if (fmpz_equal_ui(g, 1)) return;
    fmpz_gcd(g, g, &v[i]);
  }
}

// extractMatrixGCDHelper(v, g, s, A)
// 
// Input:
//   A, n x m mod s matrix.
//   n, m, s, positive integer.
// 
// Post: v, g s.t. Av = gcd(A), g = Av.
//
void extractMatrixGCDHelper(fmpz *v, fmpz *g, fmpz_t s, fmpz_mat_t A) {
  int n, m;
  n = fmpz_mat_nrows(A);
  m = fmpz_mat_ncols(A);

  fmpz_mat_t AT;
  fmpz *w;
  fmpz_t b, tmp, ci;

  fmpz_init(b);
  fmpz_init(ci);
  fmpz_init(tmp);
  fmpz_mat_init(AT, m, n);

  fmpz_mat_transpose(AT, A);
  fmpz_mat_content(b, A);
  _fmpz_vec_set(g, fmpz_mat_entry(AT, 0, 0), n);
  fmpz_set_ui(&v[0], 1);

  for (int i = 1; i < m; ++i){
    fmpz_vec_gcd(tmp, g, n);
    if (fmpz_equal(tmp, b)) break;
    w = fmpz_mat_entry(AT, i, 0);
    extractGCD(ci, g, w, s, n);
    _fmpz_vec_scalar_addmul_fmpz(g, w, n, ci);
    _fmpz_vec_scalar_mod_fmpz(g, g, n, s);
    fmpz_set(&v[i], ci);
  }

  fmpz_mat_clear(AT);
  fmpz_clear(b);
  fmpz_clear(tmp);
  fmpz_clear(ci);
}

// fmpz_mat_cpy(dst, src)
//
// Input: dst, src both matrices.
// 
// Pre: dst, src are assumed to have appropiate sizes.
//
// Post: dst is a (r1-r2) x (c2 - c1) submatrix of src that starts at (r1, c1) entry.
//
void fmpz_mat_cpy(fmpz_mat_t dst, int dr1, int dc1, int dr2, int dc2, fmpz_mat_t src, int sr1, int sc1, int sr2, int sc2) {
  fmpz_mat_t tmpSrc, tmpDst;
  fmpz_mat_window_init(tmpSrc, src, sr1, sc1, sr2, sc2);
  fmpz_mat_window_init(tmpDst, dst, dr1, dc1, dr2, dc2);

  fmpz_mat_set(tmpDst, tmpSrc);

  fmpz_mat_window_clear(tmpSrc);
  fmpz_mat_window_clear(tmpDst);
}

// fmpz_rescale(c, u, g, a, s)
//
// Input:
//   a, mod s integer.
//   s, positive integer.
//
// Post: c, u, g
// Such that (a * e) = (a, s). Where e = u + c * s/g and e is coprime to s.
//
void fmpz_rescale(fmpz_t c, fmpz_t u, fmpz_t g, fmpz_t a, fmpz_t s) {
  fmpz_t v;
  fmpz_init(v);

  fmpz_xgcd(g, u, v, a, s);
  fmpz_divexact(v, s, g);
  fmpz_stab(c, u, v, s);
  fmpz_clear(v);
}

// fmpz_rescale(ret, a, s)
//
// Input:
//   a, mod s integer.
//   s, positive integer.
//
// Post: ret = e
// Such that (a * e) = (a, s). Where e = u + c * s/g and e is coprime to s.
//
void rescaleE(fmpz_t ret, fmpz_t a, fmpz_t s) {
  fmpz_t c, u, g;

  fmpz_init(c);
  fmpz_init(u);
  fmpz_init(g);

  fmpz_rescale(c, u, g, a, s);
  fmpz_divexact(ret, s, g);
  fmpz_mul(ret, ret, c);
  fmpz_add(ret, ret, u);

  fmpz_clear(c);
  fmpz_clear(u);
  fmpz_clear(g);
}

void fmpz_mat_setRow(fmpz_mat_t dst, int row, fmpz *vec, int n) {
  for (int i = 0; i < n; ++i) {
    fmpz_set(fmpz_mat_entry(dst, row, i), &vec[i]);
  }
}

void fmpz_mat_setCol(fmpz_mat_t dst, int col, fmpz *vec, int n) {
  for (int i = 0; i < n; ++i) {
    fmpz_set(fmpz_mat_entry(dst, i, col), &vec[i]);
  }
}

// dst is assumed to be m x n
void fmpz_mat_vec_vec_mul(fmpz_mat_t dst, int m, int n, fmpz *u, fmpz *v) {
  for (int i = 0; i < m; ++i) {
    _fmpz_vec_scalar_mul_fmpz(fmpz_mat_entry(dst, i, 0), v, n, &u[i]);
  }
}

void fmpz_mat_scalar_pmod(fmpz_mat_t A, fmpz_mat_t B, fmpz_t mod) {
  int n = fmpz_mat_nrows(B);
  int m = fmpz_mat_ncols(B);
  for (int i = 0; i < n; ++i) {
    _fmpz_vec_scalar_mod_fmpz(fmpz_mat_entry(A, i, 0), fmpz_mat_entry(B, i, 0), m, mod);
  }
}

void fmpz_mat_multiplyHelper(fmpz_mat_t ret, fmpz_mat_t A, int ar1, int ac1, int ar2, int ac2, fmpz_mat_t B, int br1, int bc1, int br2, int bc2) {
  fmpz_mat_t Aw, Bw;

  fmpz_mat_window_init(Aw, A, ar1, ac1, ar2, ac2);
  fmpz_mat_window_init(Bw, B, br1, bc1, br2, bc2);

  fmpz_mat_mul(ret, Aw, Bw);

  fmpz_mat_window_clear(Aw);
  fmpz_mat_window_clear(Bw);
}

void fmpz_mod_mat_multiplyHelper(fmpz_mat_t ret, fmpz_mat_t A, int ar1, int ac1, int ar2, int ac2, fmpz_mat_t B, int br1, int bc1, int br2, int bc2, const fmpz_mod_ctx_t ctx) {
  fmpz_mat_t Aw, Bw;

  fmpz_mat_window_init(Aw, A, ar1, ac1, ar2, ac2);
  fmpz_mat_window_init(Bw, B, br1, bc1, br2, bc2);

  fmpz_mod_mat_mul(ret, Aw, Bw, ctx);

  fmpz_mat_window_clear(Aw);
  fmpz_mat_window_clear(Bw);
}


void fmpz_max(fmpz_t c, fmpz_t a, fmpz_t b) {
  if (fmpz_cmp(a, b) < 0) {
    fmpz_set(c, b);
  } else {
    fmpz_set(c, a);
  }
}
