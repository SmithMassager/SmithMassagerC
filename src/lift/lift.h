#pragma once

#include <inttypes.h>

#include "gmp.h"
#include "mpz_matrix.h"
#include "rns_matrix.h"

/** \file lift.h
 * \brief Double-plus-one lifting.
 */
struct lift_tmp {
  long * Xinv;
  rnsMatrix_t * Ap;

  rnsMatrix_t * Mp;
  rnsMatrix_t * Rx;
  rnsMatrix_t * Tp;
  rnsMatrix_t * Tx;
};
typedef struct lift_tmp lift_tmp_t;

struct lift_info_t {
  long n;
  // Number of times lift_info have been lifted. 
  int iter;
  basis_t * P;
  basis_t * X;

  lift_tmp_t tmp;

  rnsMatrix_t * Cx;
  rnsMatrix_t * Rp;
  rnsMatrix_t * Mx;
};

typedef struct lift_info_t lift_info_t;

lift_info_t * initLift(mpzMatrix_t const * A);
lift_info_t * lift(lift_info_t * info);
void finiLift(lift_info_t * info);

long numLiftIters(mpzMatrix_t const * A, mpz_t const XX);
long numLiftItersFromBound(mpz_t const XX, mpz_t const bound);

mpz_srcptr liftInfo_modulus(lift_info_t const * info);
rnsMatrix_t * liftInfo_adoptInverse(lift_info_t * info);
rnsMatrix_t * liftInfo_adoptR(lift_info_t * info);
rnsMatrix_t const * liftInfo_inverse(lift_info_t const * info);
rnsMatrix_t const * liftInfo_R(lift_info_t const * info);
rnsMatrix_t const * liftInfo_M(lift_info_t const * info);
int liftInfo_iter(lift_info_t const * info);

