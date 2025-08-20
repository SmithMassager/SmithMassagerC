#pragma once

#include "fmpz.h"
#include "fmpz_mat.h"
#include "mpz_matrix.h"

void fmpzMat_fromMpzMat(fmpz_mat_t Af, mpzMatrix_t *A);
mpzMatrix_t* fmpzMat_toMpzMat(fmpz_mat_t Af);
