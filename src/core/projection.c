#include <unistd.h>
#include <assert.h>

#include "basic.h"
#include "projection.h"

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

  //fmpz_mat_print_pretty(P);
  fmpz_mat_window_init(P1, P, 0, 0, n, r);
  fmpz_mat_window_init(E, P, 0, r, n, r+k);
  //fmpz_mat_print_pretty(E);
  //fmpz_mat_init(U, min(r, n), n);
  //fmpz_mat_init(M, n, min(r, n));
  //fmpz_mat_init(T, min(r, n), min(r, n));
  //fmpz_mat_init(F, n, k+1);
  //fmpz_mat_init(S, min(r, n), 1);

  
  //printf("%d\n", s);
  //printf("Calling SNF\n");
  for (int i = 0; i < fmpz_mat_nrows(T); ++i) {
    fmpz_set_ui(fmpz_mat_entry(T, i, i), 1);
    fmpz_set_ui(fmpz_mat_entry(S, i, 0), 1);
  }
  int success = SNF(S, U, M, T, P1, E, n, r, s, k, 0, r-1);
  //printf("Finish Calling SNF\n");
  fmpz_mat_neg(M, M);

  fmpz_mat_window_clear(P1);
  fmpz_mat_window_clear(E);
  //fmpz_mat_clear(T);
  //fmpz_mat_clear(F);
  //printf("computeProjBasis returning\n");
  return success;
}

int SNF(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t M, fmpz_mat_t T, fmpz_mat_t P, fmpz_mat_t E, int n, int r, fmpz_t s, int k, int lower, int upper) {
  //printf("SNF sarting to execute\n");
  int success = 1;
  //printf("SNF executing if statements \n");
  //printf("%d\n", s);
  //fmpz_print(s);
  if (fmpz_cmp_ui(s, 1) == 0) return 1;
  if (fmpz_cmp_ui(s, 0) == 0) return 1;
  //printf("SNF finish executing if statements \n");

  //printf("SNF initializing var\n");
  int l = fmpz_mat_nrows(S);
  int j = l - 1 - lower;
  int mid = (lower+upper)/2;
  fmpz_t tmp;
  fmpz_init(tmp);
  //printf("SNF finish initializing var\n");
  //fflush(stdout);

  if (lower == upper) {
    //printf("In base case\n");
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

    //printf("Base case: initializing var in\n");
    //printf("k: %d\n", k);
    fmpz_mat_cpy(F, 0, 0, n, 1, P, 0, lower, n, lower+1);
    fmpz_mat_cpy(F, 0, 1, n, k+1, E, 0, 0, n, k);
    //fmpz_mat_print_pretty(E);
    //fmpz_mat_print_pretty(F);
    //printf("Base case: extractMatrixGCDHelper \n");
    extractMatrixGCDHelper(v, p, s, F);
    //printf("p: ");
    //_fmpz_vec_print(p, n);
    //printf("\n");
    //fflush(stdout);
    //fprintf(stderr,"Base case: vectorGCD \n");
    vectorGCD(u, s, p, n);
    //printf("u: ");
    //_fmpz_vec_print(u, n);
    //printf("\n");
    _fmpz_vec_dot(tmp, u, p, n);
    //fprintf(stderr, "Base case: rescaleE \n");
    rescaleE(e, tmp, s);
    fmpz_mod(e, e, s);
    _fmpz_vec_scalar_mul_fmpz(p, p, n, e);
    _fmpz_vec_scalar_mod_fmpz(p, p, n, s);
    _fmpz_vec_dot(up, u, p, n);
    //printf("up: ");
    //fmpz_print(up);
    //printf("\n");
    fmpz_mod(up, up, s);
    //printf("up mod s: ");
    //fmpz_print(up);
    //printf("\n");
    if (fmpz_cmp_ui(up, 0) == 0) goto SNFclear;
    // Calculate modp(p/up, s/up) and modp(u, s/up);
    fmpz_divexact(snew, s, up);
    //fmpz_invmod(tmp, up, snew);
    //_fmpz_vec_scalar_mul_fmpz(p, p, n, tmp);
    //printf("After and before divide");
    //_fmpz_vec_print(p, n);
    _fmpz_vec_scalar_divexact_fmpz(p, p, n, up);
    //_fmpz_vec_print(p, n);
    _fmpz_vec_scalar_mod_fmpz(p, p, n, snew);
    _fmpz_vec_scalar_mod_fmpz(u, u, n, snew);
    // Set U, M to contain u, p respectively.
    fmpz_mat_setRow(U, j, u, n);
    fmpz_mat_setCol(M, j, p, n);
    //printf("setting S[j] to ");
    //fmpz_print(snew);
    //printf("\n");
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
    //printf("E:");
    //fmpz_mat_print_pretty(E);
    //printf("\n");
    //printf("T:");
    //fmpz_mat_print_pretty(T);
    //printf("\n");
    //printf("U:");
    //fmpz_mat_print_pretty(U);
    //printf("\n");
    //printf("M:");
    //fmpz_mat_print_pretty(M);
    //printf("\n");
    //printf("S:");
    //fmpz_mat_print_pretty(S);
    //printf("\n");

SNFclear:
    _fmpz_vec_clear(v, k+1);
    _fmpz_vec_clear(p, n);
    _fmpz_vec_clear(u, n);
    fmpz_clear(up);
    fmpz_clear(snew);
    fmpz_clear(e);
    fmpz_mat_clear(muE);
    fmpz_mat_clear(F);
    //printf("all cleared\n");
  } else {
    //printf("In Divide and Conquer phase\n");
    success = SNF(S, U, M, T, P, E, n, r, s, k, lower, mid);
    //printf("Success: %d\n", success);
    if (!success) goto SNFend;
    fmpz_t divisor;
    fmpz_mod_ctx_t sctx, dctx;
    fmpz_mod_mat_t Q, TW, MQUP, UP, PW, PWCpy;

    //printf("D&C: Initialize var\n");
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

    //printf("D&C: Updating P, E\n");
    //printf("D&C: Q size, %d %d\n", fmpz_mat_nrows(Q), fmpz_mat_ncols(Q));
    //printf("D&C: T size, %d %d\n", fmpz_mat_nrows(T), fmpz_mat_ncols(T));
    //printf("D&C: TW size, %d %d\n", fmpz_mat_nrows(TW), fmpz_mat_ncols(TW));
    fmpz_mod_mat_inv(Q, TW, sctx);
    //printf("D&C: Finish calculating Q\n");
    //fmpz_mat_print_pretty(Q);
    fmpz_mod_mat_multiplyHelper(UP, U, l-1-mid, 0, l-lower, n, P, 0, mid+1, n, upper+1, sctx);
    //printf("D&C: Finish calculating UP\n");
    //fmpz_mat_print_pretty(UP);
    fmpz_mod_mat_mul(UP, Q, UP, sctx);
    //fmpz_mod_mat_multiplyHelper(UP, Q, 0, 0, fmpz_mat_nrows(Q), fmpz_mat_ncols(Q), UP, 0, 0, fmpz_mat_nrows(UP), fmpz_mat_ncols(UP), sctx);
    //printf("D&C: Finish calculating QUP\n");
    //fmpz_mat_print_pretty(UP);
    fmpz_mod_mat_multiplyHelper(MQUP, M, 0, l-1-mid, n, l-lower, UP, 0, 0, fmpz_mat_nrows(UP), fmpz_mat_ncols(UP), sctx);
    //printf("D&C: Finish calculating Q, MQUP\n");
    //fmpz_mat_print_pretty(MQUP);
    //fmpz_mod_mat_get_fmpz_mat(MQUP, MQUP, sctx);
    //printf("D&C: Subtracing P - MQUP\n");
    fmpz_mod_mat_sub(PW, PW, MQUP, sctx);
    //printf("D&C: Getting divisor\n");
    //fmpz_print(fmpz_mat_entry(S, l-1-mid, 0));
    if (fmpz_cmp_ui(fmpz_mat_entry(S, l-1-mid, 0), 0) == 0) {
      fmpz_mat_zero(PW);
      goto SNFClearElse;
    }
    fmpz_divexact(divisor, s, fmpz_mat_entry(S, l-1-mid, 0));
    //printf("divisor is: ");
    //fmpz_print(divisor);
    //printf("\n");
    fmpz_mod_ctx_set_modulus(dctx, divisor);
    //printf("D&C: Setting PWCpy\n");
    //printf("D&C: PW\n");
    //fmpz_mat_print_pretty(PW);
    fmpz_mod_mat_set_fmpz_mat(PWCpy, PW, dctx);
    //printf("D&C: PW mod dctx\n");
    //fmpz_mat_print_pretty(PWCpy);
    if (!fmpz_mat_is_zero(PWCpy)) {
      printf("D&C: P - MQUP not divisble\n");
      success = 0;
      goto SNFClearElse;
    }
    fmpz_mat_scalar_divexact_fmpz(PW, PW, divisor);
    //printf("D&C: Finish Updating P, E\n");
    
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
  //printf("SNF returning\n");
  if (success && lower != upper) {
    return SNF(S, U, M, T, P, E, n, r, fmpz_mat_entry(S, l-1-mid, 0), k, mid+1, upper);
  }
  else return success;
}
