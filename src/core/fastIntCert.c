#include <stdlib.h>

#include "fastIntCert.h"
#include "rns_matrix.h"
#include "lift.h"
#include "highorder.h"
#include "imlsolve.h"
#include "reconstruct.h"
#include "fmpz_conv.h"

int fastIntCert(mpz_t s, mpzMatrix_t *Res, mpzMatrix_t *A, mpzMatrix_t *B, int n, int m) {
  rnsMatrix_t *Ap, *Bp, *Cp;
  basis_t *P;
  mpz_ptr h;
  mpz_t thres, BMX, tmp, bound, AMX;
  mpzMatrix_t *sBR, *R;
  int nprimes, k;
  long *primes;

  sBR = mpzMatrix_init(B->nrows, B->ncols);
  n = A->nrows;
  mpz_init_set_ui(thres, 1);
  mpz_init(BMX);
  mpz_init(AMX);
  mpz_init_set_ui(bound, 1);
  mpz_init(tmp);

  if (!A || !B) { return 0; }

  // Claculate thres >= 1.2 * s * n^2 * ||A|| * ||B||
  mpz_mul(thres, thres, s);
  mpzMatrix_max(AMX, A);
  mpz_mul(thres, thres, AMX);
  mpzMatrix_max(BMX, B);
  mpz_mul(thres, thres, BMX);
  mpz_mul_ui(thres, thres, 12);
  mpz_mul_ui(thres, thres, n);
  mpz_mul_ui(thres, thres, n);
  mpz_cdiv_q_ui(thres, thres, 10);

  //R = highOrderResidue(A);
  lift_info_t * info;
  // Bound = 2 s n^(n/2) ||A||^(n-1) ||B||
  mpz_pow_ui(tmp, AMX, n-1);
  mpz_ui_pow_ui(bound, n, n/2+1);
  mpz_mul(bound, bound, BMX);
  mpz_mul_ui(bound, bound, 2);
  mpz_mul(bound, bound, s);
  mpz_mul(bound, bound, tmp);
  //info = horBaseFromBound(A, NULL, bound);
  info = horBase(A, NULL);
  R = mpzMatrix_init(n, n);
  mpzMatrix_reconstruct(R, liftInfo_R(info));
  h = liftInfo_modulus(info);
  k = liftInfo_iter(info);

  // Calculate sRB
  //mpzMatrix_gemm(Res, R, B);
  mpzMatrix_rnsGemm(sBR, R, B);
  mpzMatrix_scale(sBR, s);
  int returnFalse;
  mpz_set_ui(tmp, 10); mpz_pow_ui(tmp, tmp, 60);

  if (mpz_cmp(s, tmp) <= 0) {
    primes = genCoPrimes(pickStartModulus(n), thres, &nprimes, s);
    P = basis_init(primes, nprimes);
    Ap = rnsMatrix_init(n, n, P);
    Bp = rnsMatrix_init(B->nrows, B->ncols, P);
    Cp = rnsMatrix_init(B->nrows, B->ncols, P);
    rnsMatrix_fromMpzMatrix(Ap, A);
    rnsMatrix_fromMpzMatrix(Bp, sBR);
    rnsMatrix_inverse(Ap);
    rnsMatrix_gemm(Cp, Ap, Bp);
    mpzMatrix_reconstruct(Res, Cp);
    //mpzMatrix_print(stdout, Res);

    // Check the answer A^(-1)sRB is integral iff ||rem(A^(-1)sRB, thres)||  < 0.6*s*n*||B||
    mpzMatrix_max(tmp, Res);
    mpz_mul(BMX, BMX, s);
    mpz_mul_ui(BMX, BMX, 6);
    mpz_mul_ui(BMX, BMX, n);
    mpz_cdiv_q_ui(BMX, BMX, 10);
    returnFalse = (mpz_cmp(tmp, BMX) >= 0);
    if (!returnFalse)  {
      // Calculate Res * h^k mod s
      mpzMatrix_modp(Res, s);
      mpz_powm_ui(h, h, k, s);
      mpzMatrix_scale(Res, h);
      mpzMatrix_modp(Res, s);
    }
    rnsMatrix_fini(Ap);
    rnsMatrix_fini(Bp);
    rnsMatrix_fini(Cp);
    basis_fini(P);
    free(primes);
  } else {
    imlSolve(Res, tmp, A, sBR);
    returnFalse = (mpz_cmp_ui(tmp, 1) != 0);
    if (!returnFalse) { mpzMatrix_modp(Res, s); }
  }

  mpzMatrix_fini(R);
  mpzMatrix_fini(sBR);
  mpz_clear(thres);
  mpz_clear(BMX);
  mpz_clear(AMX);
  mpz_clear(tmp);
  mpz_clear(bound);
  finiLift(info);

  return (!returnFalse);
}


int fmpzFastIntCert(fmpz_t s, fmpz_mat_t Res, fmpz_mat_t A, fmpz_mat_t B, int n, int m) {
  mpzMatrix_t *AA, *BB, *ResR;
  mpz_t ss;
  mpz_init(ss);
  fmpz_get_mpz(ss, s);
  AA = fmpzMat_toMpzMat(A);
  BB = fmpzMat_toMpzMat(B);

  ResR = mpzMatrix_init(fmpz_mat_nrows(A), fmpz_mat_ncols(B));

  int success = fastIntCert(ss, ResR, AA, BB, n, m);

  if (success) fmpzMat_fromMpzMat(Res, ResR);

  mpz_clear(ss);
  mpzMatrix_fini(AA);
  mpzMatrix_fini(BB);
  mpzMatrix_fini(ResR);

  return success;
}
