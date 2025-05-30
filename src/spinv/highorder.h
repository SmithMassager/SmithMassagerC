#pragma once

#include "rns_matrix.h"
#include "mpz_matrix.h"
#include "lift.h"

/** \file highorder.h
 * \brief High order residue computation.
 */

mpzMatrix_t * highOrderResidue_mpz(mpzMatrix_t const * A, rnsMatrix_t ** Ainv);
rnsMatrix_t * highOrderResidue_rns(mpzMatrix_t const * A, rnsMatrix_t ** Ainv);

mpzMatrix_t * highOrderResidue(mpzMatrix_t const * A);

lift_info_t * horBase(mpzMatrix_t const * A, rnsMatrix_t ** Ainv);
lift_info_t * horBaseHelper(mpzMatrix_t const * A, rnsMatrix_t ** Ainv, lift_info_t *info, int k);
lift_info_t * horBaseFromBound(mpzMatrix_t const * A, rnsMatrix_t ** Ainv, mpz_t bound);
