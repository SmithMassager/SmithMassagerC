#pragma once

#include "mpz_matrix.h"
#include "rns_matrix.h"

/** \file rns_conversion.h
 * \brief Conversion between residue number systems.
 */

#ifndef THREAD
  #define rnsMatrix_convert rnsMatrix_convertFancy
#endif

void rnsMatrix_convertFancy(rnsMatrix_t * Aq, rnsMatrix_t const * Ap);
void rnsMatrix_convertSimple(rnsMatrix_t * Aq, rnsMatrix_t const * Ap);

