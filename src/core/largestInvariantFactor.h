#pragma once

#include "fmpz_mat.h"
#include "fmpz.h"

int largestInvariantFactor(fmpz_t s, fmpz_mat_t X, fmpz_mat_t A, fmpz_t startDim);
