#include <assert.h>

#include "fmpz.h"
#include "fmpz_mat_helper.h"
#include "fmpz_conv.h"
#include "imlsolve.h"

void fmpz_mat_cmod(fmpz_mat_t A, fmpz_mat_t B, fmpz_mat_t S) {
  assert(A->c == B->c);
  assert(A->r == B->r);
  assert(A->c == S->r);

  for (int i = 0; i < A->r; ++i) {
    for (int j = 0; j < A->c; ++j) {
      fmpz_mod(fmpz_mat_entry(A, i, j), fmpz_mat_entry(B, i, j), fmpz_mat_entry(S, j, 0));
    }
  }
}

void fmpz_mat_mod(fmpz_mat_t A, fmpz_mat_t B, fmpz_t s) {
  assert(A->c == B->c);
  assert(A->r == B->r);

  for (int i = 0; i < A->r; ++i) {
    for (int j = 0; j < A->c; ++j) {
      fmpz_mod(fmpz_mat_entry(A, i, j), fmpz_mat_entry(B, i, j), s);
    }
  }
}

void fmpz_mat_mods(fmpz_mat_t A, fmpz_mat_t B, fmpz_t s) {
  assert(A->c == B->c);
  assert(A->r == B->r);

  for (int i = 0; i < A->r; ++i) {
    for (int j = 0; j < A->c; ++j) {
      fmpz_smod(fmpz_mat_entry(A, i, j), fmpz_mat_entry(B, i, j), s);
    }
  }
}

void fmpz_mat_abs(fmpz_mat_t A, fmpz_mat_t B) {
  assert(A->c == B->c);
  assert(A->r == B->r);

  for (int i = 0; i < A->r; ++i) {
    for (int j = 0; j < A->c; ++j) {
      fmpz_abs(fmpz_mat_entry(A, i, j), fmpz_mat_entry(B, i, j));
    }
  }
}

void fmpz_mat_max(fmpz_t mx, fmpz_mat_t AA) {
  fmpz_mat_t A;
  fmpz_mat_init(A, AA->r, AA->c);
  fmpz_mat_abs(A, AA);
  fmpz_set_ui(mx, 0);
  for (int i = 0; i < A->r; ++i) {
    for (int j = 0; j < A->c; ++j) {
      fmpz_max(mx, mx, fmpz_mat_entry(A, i, j));
    }
  }
  fmpz_mat_clear(A);
}

void fmpz_mat_imlSolve(fmpz_mat_t X, fmpz_t d, fmpz_mat_t A, fmpz_mat_t B) {
  mpz_t dd;
  mpzMatrix_t *XX;

  mpzMatrix_t *AA = fmpzMat_toMpzMat(A);
  mpzMatrix_t *BB = fmpzMat_toMpzMat(B);
  mpz_init(dd);
  XX = mpzMatrix_init(A->r, B->c);

  imlSolve(XX, dd, AA, BB);


  fmpzMat_fromMpzMat(X, XX);
  fmpz_set_mpz(d, dd);

  mpzMatrix_fini(AA);
  mpzMatrix_fini(BB);
  mpzMatrix_fini(XX);
  mpz_clear(dd);
}

// Extract the smallest number, a, such that a*X/d is integral.
// Which should be  a = d/gcd(X, d)
void fmpz_mat_extractLCM(fmpz_t a, fmpz_mat_t X, fmpz_t d) {
  fmpz_mat_content(a, X);
  fmpz_gcd(a, a, d);
  fmpz_divexact(a, d, a);
}


void fmpz_mat_modDiagInv(fmpz_mat_t res, fmpz_t p, fmpz_mat_t src) {
  assert(res->r == src->r);
  assert(res->c == 1);
  assert(src->c == 1);

  fmpz_t d;
  fmpz_init(d);

  for (int i = 0; i < src->r; ++i) {
    fmpz_mod(d, fmpz_mat_entry(src, i, 0), p);
    fmpz_gcdinv(d, fmpz_mat_entry(res, i, 0), d, p);
  }

  fmpz_clear(d);
}


void fmpz_mat_modDiagMul(fmpz_mat_t res, fmpz_t p, fmpz_mat_t A, fmpz_mat_t D) {
  assert(res->r == A->r);
  assert(D->r == A->c);
  assert(res->c == D->r);
  assert(1 == D->c);

  for (int i = 0; i < A->r; ++i) {
    for (int j = 0; j < A->c; ++j) {
      fmpz_mul(fmpz_mat_entry(res, i, j), fmpz_mat_entry(A, i, j), fmpz_mat_entry(D, j, 0));
    }
  }
  //for (int i = 0; i < D->r; ++i) {
  //  fmpz_mat_t resi, Ai;
  //  fmpz_mat_window_init(resi, res, 0, i, res->r, i+1);
  //  fmpz_mat_window_init(Ai, A, 0, i, A->r, i+1);

  //  fmpz_mat_scalar_mul_fmpz(resi, Ai, fmpz_mat_entry(D, i, 0));

  //  fflush(stdout);
  //  fmpz_mat_window_clear(resi);
  //  fmpz_mat_window_clear(Ai);
  //}

  fmpz_mat_mod(res, res, p);
}


void fmpz_mat_concat_vertical3(fmpz_mat_t A, fmpz_mat_t B, fmpz_mat_t C, fmpz_mat_t D) {
  assert(A->r == B->r + C->r + D->r);
  assert(A->c == B->c);
  assert(A->c == C->c);
  assert(A->c == D->c);
  fmpz_mat_t Asub;

  fmpz_mat_window_init(Asub, A, 0, 0, B->r + C->r, A->c);

  fmpz_mat_concat_vertical(Asub, B, C);
  fmpz_mat_concat_vertical(A, Asub, D);

  fmpz_mat_window_clear(Asub);
}


int fmpz_uniCert(fmpz_mat_t A) {
  mpzMatrix_t *AA = fmpzMat_toMpzMat(A);
  
  int success = uniCert(AA);

  mpzMatrix_fini(AA);
  return success;
}
