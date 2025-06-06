#include "mpz_matrix.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "reconstruct.h"
#include "rns_matrix.h"
#include "timer.h"

void mpzMatrix_extractGCD(mpzMatrix_t *v, mpzMatrix_t *g, mpzMatrix_t const * A, mpz_t const s);
void mpzMatrix_extractVectorGCD(mpz_t c, mpzMatrix_t const * v, mpzMatrix_t const * w, mpz_t const s);

mpzMatrix_t * mpzMatrix_init(long nrows, long ncols)
{
  long i;
  long nelms = nrows * ncols;
  mpzMatrix_t * M = malloc(sizeof(mpzMatrix_t));
  mpz_t * data = (mpz_t *)malloc(nelms*sizeof(mpz_t));
  for (i = 0; i < nelms; ++i) {
    mpz_init(data[i]);
  }

  M->data = data;
  M->nrows = nrows;
  M->ncols = ncols;

  return M;
}

void mpzMatrix_fini(mpzMatrix_t * M)
{
  if (M == NULL) return;
  long i;
  for (i = 0; i < mpzMatrix_numElems(M); ++i) {
    mpz_clear(M->data[i]);
  }
  free(M->data);
  M->nrows = 0;
  M->ncols = 0;
  M->data = NULL;
  free(M);
}

mpzMatrix_t const * mpzMatrix_initFromMpz(long nrows, long ncols, mpz_t const * data)
{
  mpzMatrix_t * A = malloc(sizeof(mpzMatrix_t));
  A->nrows = nrows;
  A->ncols = ncols;
  A->data = (mpz_t *)data;

  return A;
}

mpzMatrix_t * mpzMatrix_initFromMpzColumn(mpzMatrix_t const * M, long col)
{
  long nrows = M->nrows;
  mpz_t *data = malloc(sizeof(mpz_t) * nrows);
  for (int i = 0; i < nrows; ++i) mpz_set(data[i], mmget(M, i, col));

  mpzMatrix_t * A = malloc(sizeof(mpzMatrix_t));
  A->nrows = nrows;
  A->ncols = 1;
  A->data = data;

  return A;
}

// g[..][src] = A[..][col]
void mpzMatrix_copyFromMpzColumn(mpzMatrix_t *g, mpzMatrix_t const * A, long src, long dest)
{
  assert(g->nrows == A->nrows);
  assert(src < g->ncols);
  assert(dest < A->ncols);
  long nrows = g->nrows, ncols = g->ncols;
  for (int i = 0; i < nrows; ++i) {
    mpz_set(mpzMatrix_data(g)[nrows * i + dest], mpzMatrix_constData(A)[nrows * i + src]);
  }
  return;
}

mpzMatrix_t * mpzMatrix_initSet(long nrows, long ncols, long const * A_long)
{
  long i;
  mpzMatrix_t * M = malloc(sizeof(mpzMatrix_t));
  long nelms = nrows * ncols;
  mpz_t * data = (mpz_t *)malloc(nelms*sizeof(mpz_t));
  for (i = 0; i < nelms; ++i) {
    mpz_init_set_si(data[i], A_long[i]);
  }

  M->data = data;
  M->nrows = nrows;
  M->ncols = ncols;

  return M;
}

mpzMatrix_t * mpzMatrix_initFromFile(char const * filename)
{
  FILE * f = fopen(filename, "r");

  char buffer[512];
  int rc;
  char * prc;
  long i, j, nrows, ncols;
  mpz_t * data;
  mpzMatrix_t * M;

  if (!f) { return NULL; }

  prc = fgets(buffer, 512, f);
  if (!prc) { return NULL; }

  rc = fscanf(f, "%ld %ld\n", &nrows, &ncols);
  if (rc != 2) { return NULL; }

  data = (mpz_t *)malloc(nrows*ncols*sizeof(mpz_t));

  for (i = 0; i < nrows; ++i) {
    for (j = 0; j < ncols; ++j) {
      mpz_init(data[i*ncols + j]);
      rc = gmp_fscanf(f, "%Zd\n", &data[i*ncols + j]);
      if (rc != 1) { return NULL; }
    }
  }

  M = malloc(sizeof(mpzMatrix_t));
  M->data = data;
  M->nrows = nrows;
  M->ncols = ncols;
  fclose(f);

  return M;
}

void mpzMatrix_writeToFile(char const * filename, mpzMatrix_t const * M)
{
  FILE * f = fopen(filename, "w");
  long i,j;

  fprintf(f, "%%%%MatrixMarket matrix %s %s %s\n", "array", "integer", "general");
  fprintf(f, "%ld %ld\n", M->nrows, M->ncols);

  for (j = 0; j < M->ncols; ++j) {
    for (i = 0; i < M->nrows; ++i) {
      mpz_out_str(f, 10, mmgetc(M, i, j));
      fprintf(f, "\n");
    }
  }
  fclose(f);
}

mpz_t const * mpzMatrix_constData(mpzMatrix_t const * M)
{
  return (mpz_t const *)M->data;
}
mpz_t * mpzMatrix_data(mpzMatrix_t * M)
{
  return M->data;
}
mpz_ptr mmget(mpzMatrix_t * A, int row, int col)
{
  return (mpz_ptr)mmgetc(A, row, col);
}
mpz_srcptr mmgetc(mpzMatrix_t const * A, int row, int col)
{
  long i;
  assert(0 <= row);
  assert(0 <= col);
  assert(row < A->nrows);
  assert(col < A->ncols);
  i = A->ncols * row + col;
  return A->data[i];
}

long mpzMatrix_numElems(mpzMatrix_t const * M)
{
  return M->nrows * M->ncols;
}

void mpzMatrix_rand(mpzMatrix_t * M, long l)
{
  long i;
  mpz_t * data;
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, rand());

  data = mpzMatrix_data(M);
  for (i = 0; i < mpzMatrix_numElems(M); ++i) {
    mpz_urandomb(data[i], state, l);
  }
  gmp_randclear(state);
}

void mpzMatrix_set(mpzMatrix_t * dst, mpzMatrix_t const * src)
{
  long i;
  mpz_t * dst_data = mpzMatrix_data(dst);
  mpz_t const * src_data = mpzMatrix_constData(src);
  assert(dst->nrows == src->nrows && dst->ncols == src->ncols);
  for (i = 0; i < mpzMatrix_numElems(dst); ++i) {
    mpz_set(dst_data[i], src_data[i]);
  }
}

void mpzMatrix_setEntry(mpzMatrix_t *dst, long r, long c, mpz_t val) {
  long nrows = dst->nrows, ncols = dst->ncols;
  mpz_t * dst_data = mpzMatrix_data(dst);
  mpz_set(dst_data[r * ncols + c], val);
}

void mpzMatrix_setEntry_si(mpzMatrix_t *dst, long r, long c, long val) {
  long nrows = dst->nrows, ncols = dst->ncols;
  mpz_t * dst_data = mpzMatrix_data(dst);
  mpz_set_si(dst_data[r * ncols + c], val);
}


void mpzMatrix_swap(mpzMatrix_t * A, mpzMatrix_t * B)
{
  mpzMatrix_t tmp;
  tmp = *A;
  *A = *B;
  *B = tmp;
}

void mpzMatrix_zero(mpzMatrix_t * M)
{
  long i;
  mpz_t * data = mpzMatrix_data(M);
  for (i = 0; i < mpzMatrix_numElems(M); ++i) {
    mpz_set_ui(data[i], 0);
  }
}

void mpzMatrix_identity(mpzMatrix_t * M)
{
  long i;
  long n = M->nrows;
  mpz_t * data = mpzMatrix_data(M);
  assert(M->nrows == M->ncols);

  mpzMatrix_zero(M);
  for (i = 0; i < n; ++i) {
    mpz_set_ui(data[i*n + i], 1);
  }
}

int mpzMatrix_equal(mpzMatrix_t const * A, mpzMatrix_t const * B)
{
  long i;

  if(A->nrows != B->nrows) { return 0; }
  if(A->ncols != B->ncols) { return 0; }

  for (i = 0; i < mpzMatrix_numElems(A); ++i) {
    if (0 != mpz_cmp(A->data[i], B->data[i])) { return 0; }
  }

  return 1;
}

int mpzMatrix_isZero(mpzMatrix_t const * A)
{
  long i;
  mpz_t const * data = mpzMatrix_constData(A);
  for (i = 0; i < mpzMatrix_numElems(A); ++i) {
    if (mpz_cmp_si(data[i], 0) != 0) {
      return 0;
    }
  }
  return 1;
}

void mpzMatrix_max(mpz_t ret, mpzMatrix_t const * A)
{
  long i;
  mpz_t const * data = mpzMatrix_constData(A);
  mpz_set_ui(ret, 0);
  for (i = 0; i < mpzMatrix_numElems(A); ++i) {
    if (mpz_cmpabs(data[i], ret) > 0) {
      mpz_abs(ret, data[i]);
    }
  }
}

void mpzMatrix_maxColEntry(mpz_t ret, mpzMatrix_t const * A, long col)
{
  long i;
  mpz_t const * data = mpzMatrix_constData(A);
  mpz_set_ui(ret, 0);
  for (i = 0; i < A->nrows; ++i) {
    if (mpz_cmpabs(mmgetc(A, i, col), ret) > 0) {
      mpz_abs(ret, mmgetc(A, i, col));
    }
  }
}

void mpzMatrix_maxRowEntry(mpz_t ret, mpzMatrix_t const * A, long row)
{
  long i;
  mpz_t const * data = mpzMatrix_constData(A);
  mpz_set_ui(ret, 0);
  for (i = 0; i < A->ncols; ++i) {
    if (mpz_cmpabs(mmgetc(A, row, i), ret) > 0) {
      mpz_abs(ret, mmgetc(A, row, i));
    }
  }
}

void mpzMatrix_print(FILE * stream, mpzMatrix_t const * A)
{
  long i, j;
  mpz_t const * data = mpzMatrix_constData(A);
  for (i = 0; i < A->nrows; i++)
  {
    fprintf(stream, "  ");
    for (j = 0; j < A->ncols; j++) {
      mpz_out_str(stream, 10, data[i*A->ncols+j]);
      fprintf(stream, "\t");
    }
    fprintf(stream, "\n");
  }
}

void mpzMatrix_add(mpzMatrix_t * dst, mpzMatrix_t const * src)
{
  int i;
  assert(dst->nrows == src->nrows && dst->ncols == src->ncols);
  for (i = 0; i < mpzMatrix_numElems(dst); ++i) {
    mpz_add(dst->data[i], dst->data[i], src->data[i]);
  }
}

void mpzMatrix_addmul(mpzMatrix_t * dst, mpzMatrix_t const * src, mpz_t const c)
{
  int i;
  assert(dst->nrows == src->nrows && dst->ncols == src->ncols);
  for (i = 0; i < mpzMatrix_numElems(dst); ++i) {
    mpz_addmul(dst->data[i], src->data[i], c);
  }
}


void mpzMatrix_addmul_mod(mpzMatrix_t * dst, mpz_t const c, mpzMatrix_t const * src, mpz_t const mod) {
  mpzMatrix_addmul(dst, src, c);
  mpzMatrix_modp(dst, mod);
}

void calcDetBound(mpz_t bound, mpzMatrix_t const * A)
{
  long n = A->ncols;

  mpz_t S;
  mpz_init(S);

  mpzMatrix_max(S, A);
  mpz_pow_ui(S, S, n);
  mpz_ui_pow_ui(bound, n, n/2);
  mpz_mul(bound, bound, S);
  mpz_mul_ui(bound, bound, 2);
  mpz_add_ui(bound, bound, 1);

  mpz_clear(S);
}

long calcLogDetBound(mpzMatrix_t const * A) {
  long n = A->ncols;
  long exp;
  mpz_t mx;
  mpz_init(mx);

  mpzMatrix_max(mx, A);
  mpz_get_d_2exp(&exp, mx);
  mpz_clear(mx);
  return n * exp + n * log2(n/2);
}

void mpzMatrix_determinant(mpz_t det, mpzMatrix_t const * A)
{
  rnsMatrix_t *Ap;
  basis_t * P;
  mpz_t bound;
  mpz_init(bound);

  calcDetBound(bound, A);
  P = basis_initFromBound(pickStartModulus(A->ncols), bound);
  Ap = rnsMatrix_init(A->nrows, A->ncols, P);

  rnsMatrix_fromMpzMatrix(Ap, A);
  rnsMatrix_determinant(det, Ap);

  rnsMatrix_fini(Ap);
  mpz_clear(bound);
  basis_fini(P);
}


void mpzMatrix_gemm(mpzMatrix_t * dst, mpzMatrix_t const * A, mpzMatrix_t const * B)
{
  int i, j, k;
  int l, m, n;

  /* No aliased arguments. */
  assert(dst != A);
  assert(dst != B);

  /* Dimension check:
   * l x m . m x n --> l x n */
  assert(A->ncols == B->nrows);
  assert(A->nrows == dst->nrows);
  assert(B->ncols == dst->ncols);

  l = A->nrows;
  m = A->ncols;
  n = B->ncols;
  mpzMatrix_zero(dst);
  for (i = 0; i < l; ++i) {
    for (j = 0; j < n; ++j) {
      for (k = 0; k < m; ++k) {
        mpz_addmul(dst->data[i*n + j], A->data[i*m + k], B->data[k*n + j]);
      }
    }
  }
}

static void calcProductBound(mpz_t bound, mpzMatrix_t const * A, mpzMatrix_t const * B)
{
  long n = A->ncols;

  mpz_t S, T;
  mpz_inits(S, T, 0);

  mpzMatrix_max(S, A);
  mpzMatrix_max(T, B);
  mpz_mul(bound, S, T);
  mpz_mul_ui(bound, bound, 2*n);
  mpz_add_ui(bound, bound, 1);

  mpz_clears(S, T, 0);
}

void mpzMatrix_rnsGemm(mpzMatrix_t * C, mpzMatrix_t const * A, mpzMatrix_t const * B)
{
  rnsMatrix_t *Cp, *Ap, *Bp;
  basis_t * P;
  mpz_t bound;
  mpz_init(bound);

  calcProductBound(bound, A, B);
  P = basis_initFromBound(pickStartModulus(A->ncols), bound);

  Ap = rnsMatrix_init(A->nrows, A->ncols, P);
  Bp = rnsMatrix_init(B->nrows, B->ncols, P);
  Cp = rnsMatrix_init(A->nrows, B->ncols, P);

  rnsMatrix_fromMpzMatrix(Ap, A);
  rnsMatrix_fromMpzMatrix(Bp, B);
  rnsMatrix_gemm(Cp, Ap, Bp);
  mpzMatrix_reconstruct(C, Cp);

  rnsMatrix_fini(Cp);
  rnsMatrix_fini(Ap);
  rnsMatrix_fini(Bp);
  mpz_clear(bound);
  basis_fini(P);
}

static void mpz_mods_with_half(mpz_t r, mpz_t const p, mpz_t const half_p)
{
  mpz_fdiv_r(r, r, p);
  if ( 0 < mpz_cmp(r, half_p)) {
    mpz_sub(r, r, p);
  }
}

void mpzMatrix_mods(mpzMatrix_t * dst, mpz_t const p)
{
  long i;
  mpz_t half_p;
  mpz_init(half_p);
  mpz_fdiv_q_2exp(half_p, p, 1);

  for (i = 0; i < mpzMatrix_numElems(dst); ++i) {
    mpz_mods_with_half(dst->data[i], p, half_p);
  }
  mpz_clear(half_p);
}

void mpzMatrix_modp(mpzMatrix_t * dst, mpz_t const p)
{
  long i;
  for (i = 0; i < mpzMatrix_numElems(dst); ++i) {
    mpz_fdiv_r(dst->data[i], dst->data[i],p);
  }
}


void mpzMatrix_scale(mpzMatrix_t * dst, mpz_t const c)
{
  int i;
  for (i = 0; i < mpzMatrix_numElems(dst); ++i) {
    mpz_mul(dst->data[i], dst->data[i], c);
  }
}

void mpzMatrix_extractGCD(mpzMatrix_t *v, mpzMatrix_t *g, mpzMatrix_t const * A, mpz_t const s)
{
  mpz_t tmp, b;
  mpzMatrix_t *w;
  long nrows = A->nrows, ncols = A->ncols;
  assert(v->ncols == 1);
  assert(v->nrows == nrows);
  assert(g->ncols == 1);
  assert(g->nrows == nrows);
  mpzMatrix_init(nrows, 1);
  mpzMatrix_GCD(b, A, s);
  mpzMatrix_copyFromMpzColumn(g, A, 0, 0);
  mpzMatrix_setEntry_si(v, 0, 0, 1);

  for (long i = 1; i < ncols; ++i) {
    mpzMatrix_GCD(tmp, g, s);
    if (mpz_cmp(tmp, b) == 0) break;
    mpzMatrix_copyFromMpzColumn(w, A, i, 0);
    mpzMatrix_extractVectorGCD(tmp, g, w, s);
    mpzMatrix_addmul_mod(g, tmp, w, s);
  }

  mpz_clears(tmp, b, NULL);
  mpzMatrix_fini(w);
}

//void mpzMatrix_extractVectorGCD(mpz_t c, mpzMatrix_t const * v, mpzMatrix_t const * w, mpz_t const s)
//{
//  assert(v->nrows == w->nrows);
//  assert(v->nrows == 1 && w->nrows == 1);
//  long nrows = v->nrows;
//  mpz_t g, m1, m2, tmp;
//  mpz_init(g);
//  mpzMatrix_GCD(g, w, s);
//  mpz_gcd(m1, mmgetc(v, 0, 0), s);
//  mpz_gcd(m2, mmgetc(w, 0, 0), s);
//  for (long i = 1; i < nrows; ++i) {
//    if (mpz_cmp(g, m2) == 0) break;
//    stab(c, m2, mmgetc(w, i, 0), s);
//    mpz_addmul(m2, c, mmgetc(w, i, 0));
//    mpz_mod(m2, m2, s);
//    mpz_addmul(m1, c, mmgetc(v, i, 0));
//    mpz_mod(m1, m1, s);
//  }
//  stab(c, m1, m2, s);
//  mpz_clears(g, m1, m2, tmp, NULL);
//  return;
//}

void mpzMatrix_vectorGCD(mpzMatrix_t * v, mpzMatrix_t const * u, mpz_t const s)
{
  assert(v->nrows == u->nrows);
  assert(v->ncols == u->ncols);
  
  long nrows = v->nrows;
  mpz_t g, ci;
  mpz_inits(g, ci, NULL);
  mpz_t *dst_data = mpzMatrix_data(v);
  mpz_t const * src_data = mpzMatrix_constData(u);
  mpz_set(g, src_data[0]);
  mpz_mod(g, g, s);
  mpz_set_ui(dst_data[0], 1);

  for (long i = 1; i < nrows; ++i) {
    if (mpz_cmp_ui(g, 1) == 0) return;
    stab(ci, g, src_data[i], s);
    // g = g + ci * u[i] mod s
    mpz_set(dst_data[i], ci);
    mpz_addmul_mod(g, g, ci, src_data[i], s);
  }
  
  mpz_clears(g, ci, NULL);
}

//void mpzMatrix_GCD(mpz_t res, mpzMatrix_t const * A, mpz_t const s)
//{
//  mpz_set(res, s);
//  long nrows = A->nrows, ncols = A->ncols;
//  for (long i = 0; i < nrows; ++i) for (long j = 0; j < ncols; ++j) {
//    if (mpz_cmp_ui(res, 1) == 0) return;
//    mpz_gcd(res, res, mmgetc(A, i, j));
//  }
//  return;
//}

//int mpzMatrix_findMassager(mpzMatrix_t *A, mpzMatrix_t *E, mpzMatrix_t *U,
//			    mpzMatrix_t *M, mpzMatrix_t *S, mpzMatrix_t *T,
//			    mpz_t const s, long lower, long upper)
//{
//  if (mpz_cmp_ui(s, 1) == 0) return 1;
//  int ret = 1;
//
//  if (lower == upper) {
//    // Base case.
//  } else {
//    // Divide.
//    int mid = (lower+upper)/2;
//    ret = mpzMatrix_findMassager(A, E, U, M, S, T, s, lower, mid);
//    if (ret) {
//      return mpzMatrix_findMassager(A, E, U, M, S, T, mpzMatrix_constDataEntry(S, mid, 0), mid+1, upper);
//    }
//  }
//
//  return ret;
//}
