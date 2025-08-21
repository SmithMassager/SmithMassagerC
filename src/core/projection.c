#include <unistd.h>
#include <assert.h>

#include "basic.h"
#include "projection.h"
#include "timer.h"

int computeProjBasis(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t M, fmpz_mat_t T, fmpz_mat_t P, int n, int r, fmpz_t s, int k) {
  assert(S->r == r);
  assert(S->c == 1);
  assert(U->r == r);
  assert(U->c == n);
  assert(M->r == n);
  assert(M->c == r);
  assert(T->r == r);
  assert(T->c == r);

  fmpz_mat_t P1, E, F;

  fmpz_mat_window_init(P1, P, 0, 0, n, r);
  fmpz_mat_window_init(E, P, 0, r, n, r+k);

  for (int i = 0; i < fmpz_mat_nrows(T); ++i) {
    fmpz_set_ui(fmpz_mat_entry(T, i, i), 1);
    fmpz_set_ui(fmpz_mat_entry(S, i, 0), 1);
  }
  int success;
  char lbl_buffer[100];  // Buffer size: make it large enough for your label
  snprintf(lbl_buffer, sizeof(lbl_buffer), "SNF with size %d, %d %d", r-1, 0, r-1);
  static struct timeval start_real, end_real;
  static double diff_real;
  gettimeofday(&start_real, NULL);
  success = SNF(S, U, M, T, P1, E, n, r, s, k, 0, r-1);
  gettimeofday(&end_real, NULL);
  diff_real = (end_real.tv_sec - start_real.tv_sec) + (end_real.tv_usec - start_real.tv_usec) / 1000000.0;
  printf("%s: %.5f seconds (real time)\n", lbl_buffer, diff_real);

  //success = SNF(S, U, M, T, P1, E, n, r, s, k, 0, r-1);
  fmpz_mat_neg(M, M);

  fmpz_mat_window_clear(P1);
  fmpz_mat_window_clear(E);
  return success;
}

int SNF(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t M, fmpz_mat_t T, fmpz_mat_t P, fmpz_mat_t E, int n, int r, fmpz_t s, int k, int lower, int upper) {
  int success = 1;
  if (fmpz_cmp_ui(s, 1) == 0) return 1;
  if (fmpz_cmp_ui(s, 0) == 0) return 1;

  int l = fmpz_mat_nrows(S);
  int j = l - 1 - lower;
  int mid = (lower+upper)/2;
  fmpz_t tmp;
  fmpz_init(tmp);

  if (lower == upper) {
    fmpz *v, *g, *u, *uE, *p;
    fmpz_t up, snew, e;
    fmpz_mat_t muE, F;

    fmpz_init(up);
    fmpz_init(snew);
    fmpz_init(e);
    fmpz_mat_init(muE, n, k);
    fmpz_mat_init(F, n, k+1);
    uE = _fmpz_vec_init(k);
    v = _fmpz_vec_init(k+1);
    p = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);

    fmpz_mat_cpy(F, 0, 0, n, 1, P, 0, lower, n, lower+1);
    fmpz_mat_cpy(F, 0, 1, n, k+1, E, 0, 0, n, k);
    extractMatrixGCDHelper(v, p, s, F);
    vectorGCD(u, s, p, n);
    _fmpz_vec_dot(tmp, u, p, n);
    rescaleE(e, tmp, s);
    fmpz_mod(e, e, s);
    _fmpz_vec_scalar_mul_fmpz(p, p, n, e);
    _fmpz_vec_scalar_mod_fmpz(p, p, n, s);
    _fmpz_vec_dot(up, u, p, n);
    fmpz_mod(up, up, s);
    if (fmpz_cmp_ui(up, 0) == 0) goto SNFclear;
    // Calculate modp(p/up, s/up) and modp(u, s/up);
    fmpz_divexact(snew, s, up);
    _fmpz_vec_scalar_divexact_fmpz(p, p, n, up);
    _fmpz_vec_scalar_mod_fmpz(p, p, n, snew);
    _fmpz_vec_scalar_mod_fmpz(u, u, n, snew);
    // Set U, M to contain u, p respectively.
    fmpz_mat_setRow(U, j, u, n);
    fmpz_mat_setCol(M, j, p, n);
    fmpz_set(fmpz_mat_entry(S, j, 0), snew);
    // Need to calculate (E - p . u . E)/up mod s
    fmpz_mat_fmpz_vec_mul(uE, u, n, E);
    fmpz_mat_vec_vec_mul(muE, n, k, p, uE);
    fmpz_mat_sub(E, E, muE);
    fmpz_mat_scalar_pmod(E, E, s);
    fmpz_mat_scalar_divexact_fmpz(E, E, up);
    if (lower != 0) {
      fmpz_mat_t Mtmp;
      fmpz_mat_window_init(Mtmp, M, 0, j+1, fmpz_mat_nrows(M), fmpz_mat_ncols(M));
      fmpz_mod_ctx_t sctx;
      fmpz_mod_ctx_init(sctx, fmpz_mat_entry(S, fmpz_mat_nrows(S)-1, 0));
      fmpz_mod_mat_fmpz_vec_mul(fmpz_mat_entry(T, j, j+1), u, n, Mtmp, sctx);
      fmpz_mat_window_clear(Mtmp);
    }

SNFclear:
    _fmpz_vec_clear(v, k+1);
    _fmpz_vec_clear(p, n);
    _fmpz_vec_clear(u, n);
    fmpz_clear(up);
    fmpz_clear(snew);
    fmpz_clear(e);
    fmpz_mat_clear(muE);
    fmpz_mat_clear(F);
  } else {
    // TODO: DELETE THIS AFTER
    char lbl_buffer[100];  // Buffer size: make it large enough for your label
    snprintf(lbl_buffer, sizeof(lbl_buffer), "SNF with size %d, %d %d", mid-lower, lower, mid);
    static struct timeval start_real, end_real;
    static double diff_real;
    gettimeofday(&start_real, NULL);
    success = SNF(S, U, M, T, P, E, n, r, s, k, lower, mid);
    gettimeofday(&end_real, NULL);
    diff_real = (end_real.tv_sec - start_real.tv_sec) + (end_real.tv_usec - start_real.tv_usec) / 1000000.0;
    printf("%s: %.5f seconds (real time)\n", lbl_buffer, diff_real);
    // TODO: DELETE THIS AFTER

    //success = SNF(S, U, M, T, P, E, n, r, s, k, lower, mid);

    if (!success) goto SNFend;
    fmpz_t divisor;
    fmpz_mod_ctx_t sctx, dctx;
    fmpz_mod_mat_t Q, TW, MQUP, UP, PW, PWCpy;

    fmpz_init(divisor);
    fmpz_set_ui(divisor, 1);
    fmpz_mod_ctx_init(sctx, s);
    fmpz_mod_ctx_init(dctx, divisor);
    fmpz_mod_mat_init(Q, mid-lower+1, mid-lower+1, sctx);
    fmpz_mod_mat_init(MQUP, n, upper-mid, sctx);
    fmpz_mod_mat_init(UP, mid-lower+1, upper-mid, sctx);
    fmpz_mod_mat_init(PWCpy, n, upper-mid, sctx);
    fmpz_mod_mat_window_init(TW, T, l-1-mid, l-1-mid, l-lower, l-lower, sctx);
    fmpz_mod_mat_window_init(PW, P, 0, mid+1, n, upper+1, sctx);

    fmpz_mod_mat_inv(Q, TW, sctx);
    fmpz_mod_mat_multiplyHelper(UP, U, l-1-mid, 0, l-lower, n, P, 0, mid+1, n, upper+1, sctx);
    fmpz_mod_mat_mul(UP, Q, UP, sctx);
    //fmpz_mod_mat_multiplyHelper(UP, Q, 0, 0, fmpz_mat_nrows(Q), fmpz_mat_ncols(Q), UP, 0, 0, fmpz_mat_nrows(UP), fmpz_mat_ncols(UP), sctx);
    fmpz_mod_mat_multiplyHelper(MQUP, M, 0, l-1-mid, n, l-lower, UP, 0, 0, fmpz_mat_nrows(UP), fmpz_mat_ncols(UP), sctx);
    fmpz_mod_mat_sub(PW, PW, MQUP, sctx);
    if (fmpz_cmp_ui(fmpz_mat_entry(S, l-1-mid, 0), 0) == 0) {
      fmpz_mat_zero(PW);
      goto SNFClearElse;
    }
    fmpz_divexact(divisor, s, fmpz_mat_entry(S, l-1-mid, 0));
    fmpz_mod_ctx_set_modulus(dctx, divisor);
    fmpz_mod_mat_set_fmpz_mat(PWCpy, PW, dctx);
    if (!fmpz_mat_is_zero(PWCpy)) {
      printf("D&C: P - MQUP not divisble\n");
      success = 0;
      goto SNFClearElse;
    }
    fmpz_mat_scalar_divexact_fmpz(PW, PW, divisor);
    
SNFClearElse:
    fmpz_clear(divisor);
    fmpz_mod_mat_clear(Q, sctx);
    fmpz_mod_mat_window_clear(TW, sctx);
    fmpz_mod_mat_window_clear(PW, sctx);
    fmpz_mod_mat_clear(MQUP, sctx);
    fmpz_mod_mat_clear(UP, sctx);
    fmpz_mod_mat_clear(PWCpy, sctx);
    fmpz_mod_ctx_clear(sctx);
    fmpz_mod_ctx_clear(dctx);
  }

SNFend:
  fmpz_clear(tmp);
  if (success && lower != upper) {
    // TODO: DELETE THIS AFTER
    char lbl_buffer[100];  // Buffer size: make it large enough for your label
    snprintf(lbl_buffer, sizeof(lbl_buffer), "SNF with size %d, %d %d", upper-(mid+1), mid+1, upper);
    static struct timeval start_real, end_real;
    static double diff_real;
    gettimeofday(&start_real, NULL);
    success = SNF(S, U, M, T, P, E, n, r, fmpz_mat_entry(S, l-1-mid, 0), k, mid+1, upper);
    gettimeofday(&end_real, NULL);
    diff_real = (end_real.tv_sec - start_real.tv_sec) + (end_real.tv_usec - start_real.tv_usec) / 1000000.0;
    printf("%s: %.5f seconds (real time)\n", lbl_buffer, diff_real);
    return success;
    // TODO: DELETE THIS AFTER

    //return SNF(S, U, M, T, P, E, n, r, fmpz_mat_entry(S, l-1-mid, 0), k, mid+1, upper);
  }
  else return success;
}
